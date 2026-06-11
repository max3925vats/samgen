"""Two-sided (shared-S) strand fusion geometry, on a synthetic strand."""

import math
import numpy as np

from samgen.core import anchor, orient
from samgen.geometry import build_twosided_strand, generate_twosided
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
