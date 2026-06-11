import math
from samgen.core.lattice import Lattice
from samgen.design import density as D

ROUNDED = Lattice.rounded()
RHO_FULL = ROUNDED.site_density()   # ~4.639/nm^2

def test_choose_stride():
    assert D.choose_stride(RHO_FULL, 1.0) == 2     # high density -> k=2
    assert D.choose_stride(RHO_FULL, 1.0/3) == 4   # low density  -> k=4

def test_grid_options_unique_when_both_divisible():
    # box 10x10 -> 22 x 24 sites; k=2 divides both -> unique 11x12=132
    opts = D.grid_options(22, 24, 2, ROUNDED.colsep, ROUNDED.rowsep)
    assert len(opts) == 1
    o = opts[0]
    assert (o.nx, o.ny, o.count) == (11, 12, 132)
    assert math.isclose(o.box_x, 22*ROUNDED.colsep)

def test_grid_options_bracket_when_one_axis_indivisible():
    # k=4: 24%4==0 (clean), 22%4!=0 -> bracket x to 20 and 24
    opts = D.grid_options(22, 24, 4, ROUNDED.colsep, ROUNDED.rowsep)
    counts = sorted(o.count for o in opts)
    assert counts == [30, 36]                       # 5x6 and 6x6
    below = min(opts, key=lambda o: o.count)
    above = max(opts, key=lambda o: o.count)
    assert (below.ncols, above.ncols) == (20, 24)

def test_no_seam_every_option_tiles_exactly():
    # for every emitted option, the strides divide both axes -> no boundary seam
    for opts in (D.grid_options(22, 24, 2, ROUNDED.colsep, ROUNDED.rowsep),
                 D.grid_options(22, 24, 4, ROUNDED.colsep, ROUNDED.rowsep)):
        for o in opts:
            assert o.ncols % o.kx == 0 and o.nrows % o.ky == 0

def test_explicit_grid_anisotropic():
    # ligand_grid override: exactly nx x ny ligands, strides may differ per axis
    o = D.explicit_grid(22, 24, nx=8, ny=10, colsep=ROUNDED.colsep, rowsep=ROUNDED.rowsep)
    assert (o.nx, o.ny, o.count) == (8, 10, 80)
    assert o.ncols % o.kx == 0 and o.nrows % o.ky == 0   # periodic
    # kx and ky are chosen independently (anisotropic allowed)
    assert o.kx == round(22/8) and o.ky == round(24/10)


def test_density_label_after_configure():
    d = D.Density(base="coh", ligand="ch3", density=1.0)
    d.configure(kx=2, ky=2)
    # ligand on every 2nd col AND row; base elsewhere
    assert d.label(0, 0) == "ch3"
    assert d.label(0, 1) == "coh"
    assert d.label(1, 0) == "coh"
    assert d.label(2, 2) == "ch3"

def test_density_label_anisotropic():
    d = D.Density(base="coh", ligand="ch3", density=1.0)
    d.configure(kx=2, ky=3)              # different strides per axis
    assert d.label(0, 0) == "ch3"
    assert d.label(3, 2) == "ch3"        # row 3 (÷3) and col 2 (÷2)
    assert d.label(2, 2) == "coh"        # row 2 not ÷3

def test_density_label_before_configure_raises():
    d = D.Density(base="coh", ligand="ch3", density=1.0)
    import pytest
    with pytest.raises(ValueError, match="not resolved"):
        d.label(0, 0)


def test_grid_options_drops_degenerate_low_density():
    # k larger than ncols -> the 'below' multiple is 0; it must be filtered out
    opts = D.grid_options(22, 24, 24, ROUNDED.colsep, ROUNDED.rowsep)
    assert opts                                  # not empty
    assert all(o.ncols > 0 and o.nrows > 0 for o in opts)
    assert all(o.count >= 1 for o in opts)


from samgen.core.lattice import Lattice
from samgen.interactive import resolve_density_interactive
from samgen.design.density import Density
import pytest

LAT = Lattice.rounded()

def _grid0(boxx=10.0, boxy=10.0):
    return LAT.dimensions(boxx, boxy, even_cols=True)   # (22, 24)

def test_density_missing_density_raises_clear_error():
    """No 'density:' key and no 'ligand_grid:' should give a clear ValueError,
    not a cryptic TypeError from inside choose_stride."""
    d = Density("coh", "ch3")          # density=None, ligand_grid=None
    nc, nr = _grid0()
    with pytest.raises(ValueError, match="density.*ligand_grid"):
        resolve_density_interactive(d, LAT, nc, nr, is_tty=False)


def test_density_unique_no_prompt():
    d = Density("coh", "ch3", density=1.0)              # k=2 tiles -> unique
    nc, nr = _grid0()
    opt = resolve_density_interactive(d, LAT, nc, nr, is_tty=False)
    assert (opt.nx, opt.ny, opt.count) == (11, 12, 132)

def test_density_batch_requires_choice():
    d = Density("coh", "ch3", density=1.0/3)            # k=4 -> needs choice
    nc, nr = _grid0()
    with pytest.raises(ValueError, match="density_choice"):
        resolve_density_interactive(d, LAT, nc, nr, is_tty=False)

def test_density_batch_choice_above():
    d = Density("coh", "ch3", density=1.0/3, choice="above")
    nc, nr = _grid0()
    opt = resolve_density_interactive(d, LAT, nc, nr, is_tty=False)
    assert opt.count == 36 and opt.ncols == 24          # box grew in x, no seam

def test_density_interactive_prompt_pick():
    d = Density("coh", "ch3", density=1.0/3)
    nc, nr = _grid0()
    opt = resolve_density_interactive(d, LAT, nc, nr,
                                      input_fn=lambda _: "below", is_tty=True)
    assert opt.count == 30 and opt.ncols == 20

def test_density_ligand_grid_override_no_prompt():
    d = Density("coh", "ch3", density=1.0, ligand_grid=(8, 10))
    nc, nr = _grid0()
    opt = resolve_density_interactive(d, LAT, nc, nr, is_tty=False)
    assert (opt.nx, opt.ny, opt.count) == (8, 10, 80)   # explicit, no prompt


import os
from tests._synthetic import linear_strand
from samgen.geometry import generate_geometry

def test_generate_geometry_density_reproduces_132(tmp_path):
    coh = linear_strand("COH", n_chain=11)
    ch3 = linear_strand("CH3", n_chain=11)
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 10.0, "y": 10.0, "z": 10.0},
           "design": {"type": "density", "base": "base", "ligand": "ligand",
                      "density": 1.0}}
    sam = str(tmp_path / "sam.gro")
    res = generate_geometry(cfg, {"base": coh, "ligand": ch3}, out_gro=sam,
                            is_tty=False)
    assert res.manifest["counts"]["ligand"] == 132
    bx, by, _ = res.manifest["final_box"]
    assert round(bx, 3) == 10.978 and round(by, 3) == 10.368

def test_generate_geometry_density_low_grows_box(tmp_path):
    coh = linear_strand("COH", n_chain=11)
    ch3 = linear_strand("CH3", n_chain=11)
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 10.0, "y": 10.0, "z": 10.0},
           "design": {"type": "density", "base": "base", "ligand": "ligand",
                      "density": 1.0/3, "density_choice": "above"}}
    sam = str(tmp_path / "sam.gro")
    res = generate_geometry(cfg, {"base": coh, "ligand": ch3}, out_gro=sam,
                            is_tty=False)
    assert res.manifest["counts"]["ligand"] == 36
    assert round(res.manifest["final_box"][0], 3) == 11.976


def test_density_surface_passes_periodicity(tmp_path, capsys):
    coh = linear_strand("COH", n_chain=11)
    ch3 = linear_strand("CH3", n_chain=11)
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 10.0, "y": 10.0, "z": 10.0},
           "design": {"type": "density", "base": "base", "ligand": "ligand",
                      "density": 1.0}}
    res = generate_geometry(cfg, {"base": coh, "ligand": ch3},
                            out_gro=str(tmp_path / "s.gro"), is_tty=False)
    assert res.manifest["periodicity_ok"] is True


# ── Task 5: _format_density_options label ambiguity fix ──────────────────────

from samgen.interactive import _format_density_options


def test_format_density_options_four_options_unambiguous():
    """With 4 options (min, two intermediates, max), labels must be unambiguous.

    Uses ncols=21, nrows=23, k=4 which brackets both axes and produces
    counts [25, 30, 30, 36] — exactly one min and one max.
    """
    opts = D.grid_options(21, 23, 4, ROUNDED.colsep, ROUNDED.rowsep)
    assert len(opts) == 4, "precondition: exactly 4 options expected"

    formatted = _format_density_options(opts, 1.0 / 3)

    # Exactly one [below] and exactly one [above]
    assert formatted.count("[below]") == 1
    assert formatted.count("[above]") == 1

    # Intermediate options get the [ -- ] label
    assert formatted.count("[ -- ]") == 2

    # A hint to use ligand_grid for intermediate selection
    assert "ligand_grid" in formatted


def test_format_density_options_two_options_no_intermediate():
    """With 2 options (one axis indivisible), no [ -- ] labels appear."""
    # k=4, ncols=22 (divisible by 4? 22%4=2 no), nrows=24 (24%4=0 yes) -> 2 options
    opts = D.grid_options(22, 24, 4, ROUNDED.colsep, ROUNDED.rowsep)
    assert len(opts) == 2, "precondition: exactly 2 options expected"

    formatted = _format_density_options(opts, 1.0 / 3)

    assert formatted.count("[below]") == 1
    assert formatted.count("[above]") == 1
    assert "[ -- ]" not in formatted


def test_select_grid_below_above_still_returns_min_max():
    """select_grid('below'/'above') must still return the min/max-count option."""
    opts = D.grid_options(21, 23, 4, ROUNDED.colsep, ROUNDED.rowsep)
    below = D.select_grid(opts, "below")
    above = D.select_grid(opts, "above")

    min_count = min(o.count for o in opts)
    max_count = max(o.count for o in opts)

    assert below.count == min_count
    assert above.count == max_count


def test_mixed_shape_density_surface_passes_periodicity(tmp_path):
    """Regression: mixed-shape strands (different n_chain) must not cause a
    false-positive periodicity failure.

    Before the fix, per-residue centroids on a coh (n_chain=11) base + ch3
    (n_chain=6) ligand surface under a real tilt were only ~0.262 nm apart —
    below the 0.40 nm min_spacing threshold — even though all Au sites are
    0.499 nm apart.  The fix passes placement site coords to check_surface so
    the min-distance check uses the true Au-site separation.
    """
    # Base strand is longer (11 chain carbons) than the ligand (6 chain carbons)
    # so their post-tilt centroids occupy different xy positions relative to
    # the sulfur anchor — reproducing the mixed-centroid false-positive scenario.
    base = linear_strand("COH", n_chain=11)
    ligand = linear_strand("CH3", n_chain=6)
    cfg = {
        "lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
        "box": {"x": 10.0, "y": 10.0, "z": 10.0},
        "design": {"type": "density", "base": "base", "ligand": "ligand",
                   "density": 1.0},
    }
    res = generate_geometry(cfg, {"base": base, "ligand": ligand},
                            out_gro=str(tmp_path / "mixed.gro"), is_tty=False)
    assert res.manifest["periodicity_ok"] is True, (
        "Mixed-shape surface incorrectly flagged as non-periodic"
    )
