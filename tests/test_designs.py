"""Tests for Grid, MultiLigand, and generate_geometry design integration."""

import os

import pytest

from samgen.design.fraction import Fraction


def test_fraction_is_deterministic_by_seed():
    f1 = Fraction(base="b", ligand="L", fraction=0.3, seed=7)
    f2 = Fraction(base="b", ligand="L", fraction=0.3, seed=7)
    labels1 = [f1.label(r, c) for r in range(20) for c in range(20)]
    labels2 = [f2.label(r, c) for r in range(20) for c in range(20)]
    assert labels1 == labels2                       # reproducible
    frac = labels1.count("L") / len(labels1)
    assert 0.2 < frac < 0.4                          # ~target fraction


from samgen.design import make_design
from samgen.design.density import Density
from samgen.design.uniform import Uniform

def test_factory_default_is_density():
    d = make_design({"base": "b", "ligand": "L", "density": 1.0})  # no 'type'
    assert isinstance(d, Density)

def test_factory_fraction_and_uniform():
    assert isinstance(make_design({"type": "fraction", "base": "b",
                                   "ligand": "L", "fraction": 0.3}), Fraction)
    assert isinstance(make_design({"type": "uniform", "component": "x"}), Uniform)


# ── Task 1: Grid unit tests ──────────────────────────────────────────────────

from samgen.design.grid import Grid


def test_grid_label_in_range(tmp_path):
    """Grid.label() returns the mapped key for in-range (row, col) pairs."""
    pattern = tmp_path / "pattern.dat"
    pattern.write_text("0 1 0\n1 0 1\n")

    mapping = {0: "base", 1: "ligand"}
    g = Grid.from_file(str(pattern), mapping)

    # Row 0: 0 1 0
    assert g.label(0, 0) == "base"
    assert g.label(0, 1) == "ligand"
    assert g.label(0, 2) == "base"

    # Row 1: 1 0 1
    assert g.label(1, 0) == "ligand"
    assert g.label(1, 1) == "base"
    assert g.label(1, 2) == "ligand"


def test_grid_label_out_of_range_returns_base(tmp_path):
    """Out-of-range (row, col) falls back to code 0 (the base component)."""
    pattern = tmp_path / "pattern.dat"
    pattern.write_text("0 1 0\n1 0 1\n")

    g = Grid.from_file(str(pattern), {0: "base", 1: "ligand"})

    # Both row and col out of range
    assert g.label(99, 99) == "base"
    # Row in range, col out of range
    assert g.label(0, 99) == "base"
    # Row out of range only
    assert g.label(99, 0) == "base"


def test_grid_integration_with_generate_geometry(tmp_path):
    """generate_geometry with type=grid places both components and fills the box."""
    from samgen.geometry import generate_geometry
    from tests._synthetic import linear_strand

    pattern = tmp_path / "pat.dat"
    # 3 x 3 pattern with mixed 0 and 1 values
    pattern.write_text("0 1 0\n1 0 1\n0 1 0\n")

    base_mol = linear_strand("BASE", n_chain=6, with_cap=False)
    lig_mol = linear_strand("LIGD", n_chain=6, with_cap=False)

    cfg = {
        "lattice": {"rounded": True, "tilt_alpha": 0, "tilt_beta": 0},
        "box": {"x": 4.0, "y": 4.0, "z": 8.0},
        "design": {
            "type": "grid",
            "pattern": str(pattern),
            "mapping": {0: "base", 1: "ligand"},
        },
    }
    out = str(tmp_path / "surface.gro")
    result = generate_geometry(cfg, {"base": base_mol, "ligand": lig_mol},
                               out_gro=out, is_tty=False)

    counts = result.manifest["counts"]
    total = counts["base"] + counts["ligand"]
    # Every site must be assigned one component
    assert total == result.manifest["grid"]["ncols"] * result.manifest["grid"]["nrows"]
    # Both components must appear (the pattern has both 0 and 1 values)
    assert counts["base"] > 0
    assert counts["ligand"] > 0


# ── Task 2: MultiLigand unit test ────────────────────────────────────────────

from samgen.design.multilig import MultiLigand


def test_multiligand_label_maps_all_codes(tmp_path):
    """MultiLigand.label() correctly maps codes 0, 1, 2 and pads out-of-range."""
    pattern = tmp_path / "multi.dat"
    # 2 x 3 grid with codes 0, 1, 2
    pattern.write_text("0 1 2\n2 0 1\n")

    mapping = {0: "base", 1: "ligA", 2: "ligB"}
    ml = MultiLigand.from_file(str(pattern), mapping)

    # Row 0: 0 1 2
    assert ml.label(0, 0) == "base"
    assert ml.label(0, 1) == "ligA"
    assert ml.label(0, 2) == "ligB"

    # Row 1: 2 0 1
    assert ml.label(1, 0) == "ligB"
    assert ml.label(1, 1) == "base"
    assert ml.label(1, 2) == "ligA"

    # Out-of-range falls back to code 0 -> "base"
    assert ml.label(99, 99) == "base"
    assert ml.label(0, 99) == "base"
