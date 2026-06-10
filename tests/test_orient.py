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
