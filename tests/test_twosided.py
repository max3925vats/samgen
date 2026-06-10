"""Two-sided (shared-S) strand fusion geometry.

Uses the bundled ch3 example strand (configs/inputs/ch3.gro + ch3.itp): a
21-atom methyl-capped alkanethiol, so no external data is needed.
"""

import os
import math
import numpy as np

from samgen.core.molecule import Molecule
from samgen.core import anchor, orient
from samgen.geometry import build_twosided_strand, generate_twosided

ROOT = os.path.join(os.path.dirname(__file__), "..")
INPUTS = os.path.join(ROOT, "configs", "inputs")
CH3_GRO = os.path.join(INPUTS, "ch3.gro")
CH3_ITP = os.path.join(INPUTS, "ch3.itp")

# ch3: 21 atoms; methyl cap (C13 + 3 H) = 4 atoms.
#   arm A = 21 - 4 = 17;  arm B = arm A - shared S = 16;  total = 2*21 - 9 = 33.
N = 21
N_TWOSIDED = 2 * N - 9  # 33


def _ch3():
    return Molecule.from_files("CH3", CH3_GRO, CH3_ITP)


def test_atom_count_and_split(tmp_path):
    ch3 = _ch3()
    ar = anchor.autodetect_anchor(ch3)
    res = build_twosided_strand(ch3, ar.anchor_idx, ar.cap_carbon_idx,
                                str(tmp_path / "s.gro"))
    n = len(ch3.struct.atoms)  # 21
    # arm A = N - cap(4); arm B = arm A - shared S; total = 2N - 9
    assert res.n_armA == n - 4
    assert res.n_armB == n - 5
    assert res.natoms == 2 * n - 9 == N_TWOSIDED


def test_arms_antiparallel_and_c2_symmetric(tmp_path):
    ch3 = _ch3()
    ar = anchor.autodetect_anchor(ch3)
    res = build_twosided_strand(ch3, ar.anchor_idx, ar.cap_carbon_idx,
                                str(tmp_path / "s.gro"))
    c = res.strand.coords
    S = c[res.s_index]
    A, B = c[:res.n_armA], c[res.n_armA:]
    hA = A[np.argmax(np.linalg.norm(A - S, axis=1))] - S
    hB = B[np.argmax(np.linalg.norm(B - S, axis=1))] - S
    hA /= np.linalg.norm(hA); hB /= np.linalg.norm(hB)
    assert hA @ hB < -0.999  # chain axes antiparallel

    # arm B is exactly the proper 180-deg-about-x image of arm A (no S) ->
    # chirality preserved (a reflection would NOT satisfy this).
    Rx = orient.rotation_matrix("x", 180)
    armA_noS = np.delete(A, res.s_index, axis=0)
    assert np.abs(B - armA_noS @ Rx.T).max() < 1e-6


def test_both_arms_get_correct_tilt(tmp_path):
    ch3 = _ch3()
    ar = anchor.autodetect_anchor(ch3)
    res = build_twosided_strand(ch3, ar.anchor_idx, ar.cap_carbon_idx,
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
    ch3 = _ch3()
    ar = anchor.autodetect_anchor(ch3)
    cfg = {"lattice": {"a": 0.288, "tilt_alpha": 28, "tilt_beta": 53},
           "box": {"x": 4.0, "y": 4.0, "z": 12.0}}
    res, geom = generate_twosided(ch3, ar.anchor_idx, ar.cap_carbon_idx, cfg,
                                  str(tmp_path / "strand.gro"),
                                  str(tmp_path / "surface.gro"))
    strands = list(geom.manifest["counts"].values())[0]
    assert geom.manifest["natoms"] == strands * N_TWOSIDED
    assert os.path.exists(str(tmp_path / "strand.gro"))
