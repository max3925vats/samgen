import numpy as np
from samgen.core.gro import GroStructure, GroAtom
from samgen.core.lattice import Lattice
from samgen.core import periodicity as P

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
