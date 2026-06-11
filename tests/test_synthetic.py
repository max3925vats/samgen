import numpy as np
from tests._synthetic import linear_strand

def test_linear_strand_shape_and_anchor():
    mol = linear_strand("LIG", n_chain=11, with_cap=True)
    # atoms: S + n_chain carbons + cap carbon (+3 H) ; masses/bonds present
    assert mol.masses is not None and mol.bonds is not None
    assert mol.struct.atoms[mol.anchor_index].atomname.startswith("S")
    # chain lies along +z (canonical)
    z = mol.coords[:, 2]
    assert z[mol.backbone_index(9)] > z[mol.anchor_index]

def test_linear_strand_arbitrary_size():
    for m in (7, 33):
        mol = linear_strand("X", n_chain=m, with_cap=False)
        assert len(mol.struct.atoms) == m + 1   # carbons + S, no cap
