"""Topology parsing, merging, assembly, and the grompp gate.

Uses the bundled ch3/coh example strands under configs/inputs/ (pre-oriented
along z), so these run without any external data.
"""

import os
import pytest

from samgen.core.molecule import Molecule
from samgen.core.topfile import parse_top
from samgen.geometry import generate_geometry
from samgen.topology import assemble_topology, merge_topfiles
from samgen import gmx

ROOT = os.path.join(os.path.dirname(__file__), "..")
INPUTS = os.path.join(ROOT, "configs", "inputs")
CH3_ITP = os.path.join(INPUTS, "ch3.itp")
COH_ITP = os.path.join(INPUTS, "coh.itp")
CH3_GRO = os.path.join(INPUTS, "ch3.gro")
COH_GRO = os.path.join(INPUTS, "coh.gro")


def test_parse_extracts_moleculetype_and_atomtypes():
    tf = parse_top(CH3_ITP)
    assert "CH3" in tf.molecules
    assert len(tf.atomtypes) > 0
    # every atom type used by CH3 is defined in the same self-contained .itp
    assert all(t in tf.atomtypes for t in tf.atomtype_names_used("CH3"))


def test_merge_pools_components_and_unions_atomtypes():
    merged, warns = merge_topfiles({"COH": COH_ITP, "CH3": CH3_ITP})
    assert "COH" in merged.molecules and "CH3" in merged.molecules
    # union covers both: hc is ch3-only, oh/ho are coh-only
    assert {"hc", "oh", "ho"} <= set(merged.atomtypes)


def test_assemble_writes_complete_top(tmp_path):
    coh = Molecule.from_files("COH", COH_GRO, COH_ITP)
    ch3 = Molecule.from_files("CH3", CH3_GRO, CH3_ITP)
    cfg = {"lattice": {"a": 0.288, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 4.0, "y": 4.0, "z": 10.0},
           "design": {"type": "density", "base": "base", "ligand": "ligand",
                      "fraction": 0.3, "seed": 1},
           "output": {"order": ["COH", "CH3"]}}
    sam = str(tmp_path / "sam.gro")
    generate_geometry(cfg, {"base": coh, "ligand": ch3}, out_gro=sam)
    top = str(tmp_path / "topol.top")
    counts = assemble_topology(sam, itp_map={"COH": COH_ITP, "CH3": CH3_ITP},
                               order=["COH", "CH3"], out_top=top,
                               out_gro=str(tmp_path / "ro.gro"))
    text = open(top).read()
    for sec in ("[ defaults ]", "[ atomtypes ]", "[ moleculetype ]",
                "[ system ]", "[ molecules ]"):
        assert sec in text
    assert counts["COH"] + counts["CH3"] > 0


def test_missing_atomtype_raises(tmp_path):
    # an .itp whose molecule references an undefined atom type must be rejected
    itp = tmp_path / "bad.itp"
    itp.write_text(
        "[ moleculetype ]\nBAD 3\n[ atoms ]\n"
        "1 ZZ 1 BAD X1 1 0.0 12.0\n"
    )
    gro = tmp_path / "s.gro"
    gro.write_text("t\n1\n    1BAD     X1    1   0.000   0.000   0.000\n1 1 1\n")
    with pytest.raises(ValueError, match="atom types not in"):
        assemble_topology(str(gro), itp_map={"BAD": str(itp)}, order=["BAD"],
                          out_top=str(tmp_path / "t.top"))


@pytest.mark.skipif(not gmx.available(), reason="GROMACS not on PATH")
def test_grompp_accepts_assembled_topology(tmp_path):
    coh = Molecule.from_files("COH", COH_GRO, COH_ITP)
    ch3 = Molecule.from_files("CH3", CH3_GRO, CH3_ITP)
    cfg = {"lattice": {"a": 0.288, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 4.0, "y": 4.0, "z": 10.0},
           "design": {"type": "grid", "pattern": None, "mapping": {0: "base", 1: "ligand"}},
           "output": {"order": ["COH", "CH3"]}}
    # all-zero pattern -> a uniform coh base, written on the fly
    pat = tmp_path / "p.dat"
    pat.write_text("\n".join(" ".join("0" for _ in range(20)) for _ in range(20)))
    cfg["design"]["pattern"] = str(pat)
    sam = str(tmp_path / "sam.gro")
    generate_geometry(cfg, {"base": coh, "ligand": ch3}, out_gro=sam)
    # validate=True runs grompp and raises if the topology is rejected
    assemble_topology(sam, itp_map={"COH": COH_ITP, "CH3": CH3_ITP},
                      order=["COH", "CH3"], out_top=str(tmp_path / "t.top"),
                      out_gro=str(tmp_path / "ro.gro"), validate=True)
