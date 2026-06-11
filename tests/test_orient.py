# tests/test_orient.py
import warnings
from tests._synthetic import linear_strand     # Task 11
from samgen.core.anchor import backbone_head


def test_backbone_head_picks_nth_chain_carbon():
    # S(idx anchor) - C1 - C2 - ... ; cap methyl on S is excluded
    mol = linear_strand("LIG", n_chain=11, with_cap=True)
    a = mol.anchor_index            # synthetic builder exposes this
    head = backbone_head(mol, a, n_carbons=9)
    # 9th backbone carbon from S (0-based positions defined by the builder)
    assert head == mol.backbone_index(9)


def test_backbone_head_warns_on_short_chain():
    mol = linear_strand("LIG", n_chain=4, with_cap=True)
    a = mol.anchor_index
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        head = backbone_head(mol, a, n_carbons=9)
        assert any("chain" in str(x.message).lower() for x in w)
    assert head == mol.backbone_index(4)   # stops at the last available carbon


# ── Task 10 ──────────────────────────────────────────────────────────────────

import numpy as np
from samgen.core import orient
from samgen.geometry import generate_geometry


def test_density_orientation_uses_backbone_not_headgroup(tmp_path):
    # strand whose headgroup juts sideways; whole-molecule axis would be wrong,
    # backbone axis (S->C9) is along the chain.
    mol = linear_strand("LIG", n_chain=11, with_cap=True,
                         headgroup_offset=(3.0, 0.0, 0.0))   # big sideways group
    mol_raw = mol.unoriented_copy()                          # lay chain along x
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 0, "tilt_beta": 0},
           "box": {"x": 3.0, "y": 3.0, "z": 8.0},
           "design": {"type": "uniform", "component": "c"},
           "components_meta": {"c": {"anchor": mol.anchor_index + 1,
                                     "canonicalize": True, "backbone_carbons": 9}}}
    sam = str(tmp_path / "s.gro")
    res = generate_geometry(cfg, {"c": mol_raw}, out_gro=sam, is_tty=False)
    assert res.manifest["natoms"] > 0

    # Direct check: read back the first strand and verify S->C9 points +z.
    # tilt_alpha=0 tilt_beta=0, so after canonicalization the backbone is
    # exactly along z (no further rotation applied). A wrong axis (e.g. along x)
    # would give dz ≈ 0, making this assertion fail clearly.
    from samgen.core.gro import read_gro
    struct = read_gro(sam)
    res1 = [a for a in struct.atoms if a.resid == 1]
    s_atom  = next(a for a in res1 if a.atomname == "S1")
    c9_atom = next(a for a in res1 if a.atomname == "C9")
    dz = c9_atom.z - s_atom.z
    # backbone length S->C9 is 9 * 0.15 nm = 1.35 nm; after zero-tilt
    # canonicalization dz must be positive and >= 0.9 * 1.35 = 1.215 nm
    assert dz > 1.215, (
        f"backbone not along +z after canonicalization: S.z={s_atom.z:.3f}, "
        f"C9.z={c9_atom.z:.3f}, dz={dz:.3f} nm"
    )
