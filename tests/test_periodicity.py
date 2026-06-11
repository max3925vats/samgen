import numpy as np
from samgen.core.gro import GroStructure, GroAtom
from samgen.core.lattice import Lattice
from samgen.core import periodicity as P


def _mixed_surface_2site(cs: float, box: tuple) -> GroStructure:
    """Two-residue surface where each residue has TWO atoms that pull the
    centroid toward the neighbour — simulating a mixed coh/ch3 surface.

    Site 0 at x=0: anchor at (0, 0), chain atom at (0.35, 0) -> centroid (0.175, 0)
    Site 1 at x=cs: anchor at (cs, 0), chain atom at (cs-0.35, 0) -> centroid (cs-0.175, 0)
    Centroid distance ≈ 0.149 nm << 0.40 nm min_spacing (false positive without fix).
    Actual site distance = cs = 0.499 nm (valid).
    """
    atoms = [
        GroAtom(1, "BAS", "S1", 1, 0.0, 0.0, 0.0),
        GroAtom(1, "BAS", "C1", 2, 0.35, 0.0, 0.0),
        GroAtom(2, "BAS", "S1", 3, cs, 0.0, 0.0),
        GroAtom(2, "BAS", "C1", 4, cs - 0.35, 0.0, 0.0),
    ]
    return GroStructure("mixed_2site", atoms, box)

LAT = Lattice.rounded()

def _surface(strands_xy, box, ligand_flags=None):
    """One dummy atom per strand at its (x, y); resname LIG for ligands else BAS."""
    atoms = []
    for i, (x, y) in enumerate(strands_xy):
        lig = ligand_flags[i] if ligand_flags else True
        atoms.append(GroAtom(i + 1, "LIG" if lig else "BAS", "S1", i + 1, x, y, 0.0))
    return GroStructure("t", atoms, box)

def test_clean_density_grid_passes():
    # 3 x 3 ligand grid at stride-2 spacing, box = exact multiple -> uniform
    cs, rs = LAT.colsep, LAT.rowsep
    xy = [(c*2*cs, r*2*rs) for r in range(3) for c in range(3)]
    box = (3*2*cs, 3*2*rs, 5.0)
    rep = P.check_surface(_surface(xy, box), LAT, ligand_resname="LIG")
    assert rep.ok, rep.issues

def test_seam_grid_flags():
    # same grid but the box is too narrow in x -> last column is compressed at wrap
    cs, rs = LAT.colsep, LAT.rowsep
    xy = [(c*2*cs, r*2*rs) for r in range(3) for c in range(3)]
    box = (3*2*cs - cs, 3*2*rs, 5.0)        # x shortened by one colsep -> seam
    rep = P.check_surface(_surface(xy, box), LAT, ligand_resname="LIG",
                          min_spacing=0.40)
    assert not rep.ok
    assert any("seam" in m or "min strand" in m for m in rep.issues)

def test_box_not_lattice_multiple_flags():
    box = (10.0, 10.0, 5.0)                  # 10 / 0.499 is not integer
    rep = P.check_surface(_surface([(0.0, 0.0)], box), LAT)
    assert not rep.ok
    assert any("multiple" in m for m in rep.issues)


# ---------------------------------------------------------------------------
# site_xy parameter tests
# ---------------------------------------------------------------------------

def test_site_xy_uniform_grid_passes_even_with_close_centroids():
    """site_xy overrides centroid-based distance check.

    Build a surface whose residue xy-centroids are only ~0.149 nm apart
    (mimicking a mixed coh/ch3 surface where different-shaped strands have
    centroids pulled toward adjacent sites) but whose placement sites are
    0.499 nm apart (a valid Au lattice).  Without site_xy this triggers a
    false-positive warning; with site_xy it should pass.
    """
    cs = LAT.colsep   # 0.499 nm
    box = (2 * cs, LAT.rowsep * 2, 5.0)
    struct = _mixed_surface_2site(cs, box)

    # Without site_xy the centroid distance (~0.149 nm) triggers a false positive.
    rep_no_site = P.check_surface(struct, LAT, min_spacing=0.40)
    assert not rep_no_site.ok, (
        "Expected a false-positive without site_xy (centroid distances too close)")

    # With the correct site positions the min-distance check uses 0.499 nm -> OK.
    site_xy = np.array([[0.0, 0.0], [cs, 0.0]])
    rep = P.check_surface(struct, LAT, min_spacing=0.40, site_xy=site_xy)
    assert rep.ok, rep.issues


def test_site_xy_seam_still_flags():
    """A seam in site_xy (too-short box causing wrap-around clash) is still caught."""
    cs = LAT.colsep
    # 3 valid sites but the box is shrunk so the periodic image of site 0
    # is only 0.3 nm from site 2 — a genuine seam.
    site_xy = np.array([[0.0, 0.0], [cs, 0.0], [2 * cs, 0.0]])
    # Shrink x by 0.25*cs so the wrap distance from site[2] to site[0] is
    # cs - 0.25*cs = 0.75*cs ≈ 0.374 nm < 0.40 nm minimum.
    box = (3 * cs - 0.25 * cs, LAT.rowsep * 2, 5.0)
    struct = _surface(list(site_xy), box)

    rep = P.check_surface(struct, LAT, min_spacing=0.40, site_xy=site_xy)
    assert not rep.ok
    assert any("min strand" in m or "seam" in m for m in rep.issues)
