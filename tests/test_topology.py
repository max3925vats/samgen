"""Topology parsing, merging, assembly, and the grompp gate (synthetic inputs)."""

import os
import pytest

from samgen.core.molecule import Molecule
from samgen.core.topfile import parse_top
from samgen.geometry import generate_geometry
from samgen.topology import assemble_topology, merge_topfiles
from samgen import gmx
from tests._synthetic import linear_strand

# A minimal self-contained .itp covering the two synthetic resnames.
_ITP = """[ atomtypes ]
 C  C 0.0 0.0 A 3.4e-01 4.6e-01
 S  S 0.0 0.0 A 3.6e-01 1.0e+00
[ moleculetype ]
{name} 3
[ atoms ]
{atoms}
"""

def _write_itp(path, name, n_atoms):
    rows = []
    for i in range(n_atoms):
        typ = "S" if i == 0 else "C"
        mass = 32.06 if i == 0 else 14.027
        rows.append(f" {i+1} {typ} 1 {name} A{i+1} {i+1} 0.0 {mass}")
    path.write_text(_ITP.format(name=name, atoms="\n".join(rows)))
    return str(path)


def test_parse_and_merge(tmp_path):
    a = _write_itp(tmp_path / "a.itp", "AAA", 5)
    b = _write_itp(tmp_path / "b.itp", "BBB", 6)
    merged, warns = merge_topfiles({"AAA": a, "BBB": b})
    assert "AAA" in merged.molecules and "BBB" in merged.molecules
    assert {"C", "S"} <= set(merged.atomtypes)


def test_assemble_writes_complete_top(tmp_path):
    coh = linear_strand("AAA", n_chain=4, with_cap=False)   # 5 atoms
    lig = linear_strand("BBB", n_chain=5, with_cap=False)   # 6 atoms
    a = _write_itp(tmp_path / "a.itp", "AAA", 5)
    b = _write_itp(tmp_path / "b.itp", "BBB", 6)
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 0, "tilt_beta": 0},
           "box": {"x": 4.0, "y": 4.0, "z": 8.0},
           "design": {"type": "fraction", "base": "base", "ligand": "ligand",
                      "fraction": 0.3, "seed": 1},
           "output": {"order": ["AAA", "BBB"]}}
    sam = str(tmp_path / "sam.gro")
    generate_geometry(cfg, {"base": coh, "ligand": lig}, out_gro=sam, is_tty=False)
    top = str(tmp_path / "topol.top")
    counts = assemble_topology(sam, itp_map={"AAA": a, "BBB": b},
                               order=["AAA", "BBB"], out_top=top,
                               out_gro=str(tmp_path / "ro.gro"))
    text = open(top).read()
    for sec in ("[ defaults ]", "[ atomtypes ]", "[ moleculetype ]",
                "[ system ]", "[ molecules ]"):
        assert sec in text
    assert counts["AAA"] + counts["BBB"] > 0


def test_missing_atomtype_raises(tmp_path):
    itp = tmp_path / "bad.itp"
    itp.write_text("[ moleculetype ]\nBAD 3\n[ atoms ]\n1 ZZ 1 BAD X1 1 0.0 12.0\n")
    gro = tmp_path / "s.gro"
    gro.write_text("t\n1\n    1BAD     X1    1   0.000   0.000   0.000\n1 1 1\n")
    with pytest.raises(ValueError, match="atom types not in"):
        assemble_topology(str(gro), itp_map={"BAD": str(itp)}, order=["BAD"],
                          out_top=str(tmp_path / "t.top"))


# ── Task 4: grompp validate-gate with real shipped GAFF inputs ────────────────

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")


@pytest.mark.skipif(not gmx.available(), reason="gmx not on PATH")
def test_grompp_accepts_ch3_surface(tmp_path):
    """Assemble a CH3 surface from the shipped GAFF .gro/.itp and validate with grompp.

    Uses configs/inputs/ch3.gro (resname CH3, anchor S12) — a complete GAFF
    topology that grompp can fully validate without synthetic stand-ins.
    """
    from samgen.core.molecule import Molecule

    ch3_gro = os.path.join(_REPO_ROOT, "configs", "inputs", "ch3.gro")
    ch3_itp = os.path.join(_REPO_ROOT, "configs", "inputs", "ch3.itp")

    mol = Molecule.from_files("CH3", ch3_gro, ch3_itp)

    cfg = {
        "lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
        "box": {"x": 4.0, "y": 4.0, "z": 8.0},
        "design": {"type": "uniform", "component": "ch3"},
        # Canonicalize from the pre-oriented ch3.gro (anchor = S12 by name)
        "components_meta": {
            "ch3": {
                "canonicalize": True,
                "anchor": "S12",
                "backbone_carbons": 9,
            }
        },
    }

    surface_gro = str(tmp_path / "surface.gro")
    generate_geometry(cfg, {"ch3": mol}, out_gro=surface_gro, is_tty=False)

    out_top = str(tmp_path / "topol.top")
    out_gro = str(tmp_path / "reordered.gro")

    # validate=True runs gmx grompp as a hard gate; raises RuntimeError on failure
    counts = assemble_topology(
        surface_gro,
        itp_map={"CH3": ch3_itp},
        order=["CH3"],
        out_top=out_top,
        out_gro=out_gro,
        validate=True,
    )

    # If we reach here grompp accepted the topology — assert count is sane
    assert counts["CH3"] > 0
