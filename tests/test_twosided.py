"""Two-sided (shared-S) strand fusion geometry, on a synthetic strand."""

import math
import numpy as np

from samgen.core import anchor, orient
from samgen.geometry import build_twosided_strand, generate_twosided
from samgen.topology import assemble_topology
from tests._synthetic import linear_strand


def _strand():
    # S + 11 carbons + methyl cap (Cc + 3H) = 16 atoms
    return linear_strand("LIG", n_chain=11, with_cap=True)


def test_atom_count_and_split(tmp_path):
    mol = _strand()
    ar = anchor.autodetect_anchor(mol)
    res = build_twosided_strand(mol, ar.anchor_idx, ar.cap_carbon_idx,
                                str(tmp_path / "s.gro"))
    n = len(mol.struct.atoms)
    assert res.n_armA == n - 4            # cap carbon + 3 H removed
    assert res.n_armB == n - 5            # arm B drops the shared S too
    assert res.natoms == 2 * n - 9


def test_arms_antiparallel_and_c2_symmetric(tmp_path):
    mol = _strand()
    ar = anchor.autodetect_anchor(mol)
    res = build_twosided_strand(mol, ar.anchor_idx, ar.cap_carbon_idx,
                                str(tmp_path / "s.gro"))
    c = res.strand.coords
    S = c[res.s_index]
    A, B = c[:res.n_armA], c[res.n_armA:]
    hA = A[np.argmax(np.linalg.norm(A - S, axis=1))] - S
    hB = B[np.argmax(np.linalg.norm(B - S, axis=1))] - S
    hA /= np.linalg.norm(hA); hB /= np.linalg.norm(hB)
    assert hA @ hB < -0.999
    Rx = orient.rotation_matrix("x", 180)
    armA_noS = np.delete(A, res.s_index, axis=0)
    assert np.abs(B - armA_noS @ Rx.T).max() < 1e-6


def test_both_arms_get_correct_tilt(tmp_path):
    mol = _strand()
    ar = anchor.autodetect_anchor(mol)
    res = build_twosided_strand(mol, ar.anchor_idx, ar.cap_carbon_idx,
                                str(tmp_path / "s.gro"))
    t = orient.apply_tilt(res.strand.coords, alpha=28, beta=53)
    S = t[res.s_index]
    A, B = t[:res.n_armA], t[res.n_armA:]
    hA = A[np.argmax(np.linalg.norm(A - S, axis=1))] - S
    hB = B[np.argmax(np.linalg.norm(B - S, axis=1))] - S
    hA /= np.linalg.norm(hA); hB /= np.linalg.norm(hB)
    aA = math.degrees(math.acos(np.clip(hA @ [0, 0, 1.0], -1, 1)))
    aB = math.degrees(math.acos(np.clip(hB @ [0, 0, -1.0], -1, 1)))
    assert abs(aA - 28) < 0.5 and abs(aB - 28) < 0.5


def test_full_twosided_surface(tmp_path):
    """generate_twosided builds a tiled two-sided SAM surface with correct atom count."""
    mol = _strand()
    ar = anchor.autodetect_anchor(mol)

    config = {
        "lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
        "box": {"x": 4.0, "y": 4.0, "z": 12.0},
    }

    strand_gro = str(tmp_path / "strand.gro")
    surface_gro = str(tmp_path / "surface.gro")

    res, geom = generate_twosided(
        mol, ar.anchor_idx, ar.cap_carbon_idx,
        config, strand_gro, surface_gro,
    )

    # strand .gro must be written
    import os
    assert os.path.isfile(strand_gro), "strand .gro was not written"

    # total atoms = strands * atoms_per_fused_strand
    # atoms_per_fused_strand = 2 * n_one_sided - 9  (cap C + 3H + shared S saved once)
    n_one_sided = len(mol.struct.atoms)          # 16 for the synthetic strand
    atoms_per_strand = 2 * n_one_sided - 9       # == res.natoms
    strands = geom.manifest["grid"]["ncols"] * geom.manifest["grid"]["nrows"]
    assert geom.manifest["natoms"] == strands * atoms_per_strand


# ── Task 3: Integrated topology assembly from a two-sided surface ─────────────

# Minimal self-contained .itp template for the fused two-sided strand.
# Atom 0 is type S (mass 32.06); the remaining 22 are type C (mass 14.027).
_TWOSIDED_ITP_HEADER = """\
[ atomtypes ]
 C  C 0.0 0.0 A 3.4e-01 4.6e-01
 S  S 0.0 0.0 A 3.6e-01 1.0e+00
[ moleculetype ]
{name} 3
[ atoms ]
{atoms}
"""


def _write_twosided_itp(path, name: str, n_atoms: int) -> str:
    """Write a minimal .itp for a fused two-sided strand (atom 0 = S, rest C)."""
    rows = []
    for i in range(n_atoms):
        typ = "S" if i == 0 else "C"
        mass = 32.06 if i == 0 else 14.027
        rows.append(f" {i+1} {typ} 1 {name} A{i+1} {i+1} 0.0 {mass}")
    path.write_text(_TWOSIDED_ITP_HEADER.format(name=name, atoms="\n".join(rows)))
    return str(path)


def test_twosided_topology_assembly(tmp_path):
    """Build a two-sided SAM surface, assemble a merged topology, and verify it."""
    mol = linear_strand("LIG", n_chain=11, with_cap=True)
    ar = anchor.autodetect_anchor(mol)

    cfg = {
        "lattice": {"rounded": True, "tilt_alpha": 28, "tilt_beta": 53},
        "box": {"x": 4.0, "y": 4.0, "z": 12.0},
    }

    out_strand = str(tmp_path / "strand.gro")
    out_surface = str(tmp_path / "surface.gro")

    res, geom = generate_twosided(
        mol, ar.anchor_idx, ar.cap_carbon_idx,
        cfg, out_strand, out_surface,
    )

    # Fused strand must have the expected atom count (2*16 - 9 = 23)
    expected_natoms_per_strand = 2 * len(mol.struct.atoms) - 9
    assert res.natoms == expected_natoms_per_strand

    resname = res.strand.name   # "LIG"

    # Write a minimal .itp for the fused strand
    itp_path = _write_twosided_itp(
        tmp_path / f"{resname}.itp", resname, res.natoms
    )

    out_top = str(tmp_path / "topol.top")
    out_gro = str(tmp_path / "reordered.gro")

    counts = assemble_topology(
        out_surface,
        itp_map={resname: itp_path},
        order=[resname],
        out_top=out_top,
        out_gro=out_gro,
    )

    # The written topology must contain required GROMACS sections
    top_text = open(out_top).read()
    assert "[ moleculetype ]" in top_text
    assert "[ system ]" in top_text
    assert "[ molecules ]" in top_text

    # Per-resname count must match the number of placed two-sided strands
    placed = geom.manifest["counts"]["twosided"]
    assert counts[resname] == placed
