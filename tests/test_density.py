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
