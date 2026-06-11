"""Strand orientation: apply the SAM tilt, validate, optionally canonicalize.

Orientation is layered, never silent magic:
    1. Default trusts the input (chain ~ along z) and only applies the tilt.
    2. A sanity check always runs and fails LOUD if a strand is mis-oriented.
    3. Chemistry-anchored canonicalization is opt-in (anchor -> head vector),
       avoiding the sign/degeneracy failure modes of principal-axis methods.
"""

from __future__ import annotations

import math
import numpy as np


def rotation_matrix(axis: str, degrees: float) -> np.ndarray:
    """Right-handed rotation about a Cartesian axis."""
    t = math.radians(degrees)
    c, s = math.cos(t), math.sin(t)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    raise ValueError(f"axis must be x|y|z, got {axis!r}")


def apply_tilt(coords: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """Apply the SAM tilt to a canonically-oriented strand (chain along z).

    Two moves: twist `beta` about z, then tilt `alpha` about y. Rotations are
    about the centroid so we don't translate the strand. alpha = tilt from
    surface normal; beta = twist about chain axis.
    """
    centroid = coords.mean(axis=0)
    centred = coords - centroid
    rot = rotation_matrix("y", -alpha) @ rotation_matrix("z", -beta)
    return (centred @ rot.T) + centroid


def long_axis(coords: np.ndarray) -> np.ndarray:
    """Unit vector of the molecule's longest extent (PCA, for validation only)."""
    centred = coords - coords.mean(axis=0)
    # Largest-variance direction == smallest moment of inertia == chain axis.
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    return vt[0] / np.linalg.norm(vt[0])


def check_oriented(coords: np.ndarray, tol_deg: float = 20.0) -> float:
    """Return the chain-axis angle from z (deg). Raise if it exceeds tol_deg.

    This is the loud failure that prevents building a garbage surface from a
    mis-oriented input. Caller decides whether to canonicalize instead.
    """
    axis = long_axis(coords)
    cos = abs(np.dot(axis, np.array([0.0, 0.0, 1.0])))
    angle = math.degrees(math.acos(min(1.0, cos)))
    if angle > tol_deg:
        raise ValueError(
            f"strand chain axis is {angle:.1f} deg off z (tol {tol_deg}). "
            "Pre-orient the file, or enable chemistry-anchored canonicalization "
            "by specifying anchor/head atoms."
        )
    return angle


def canonicalize(coords: np.ndarray, anchor_idx: int, head_idx: int) -> np.ndarray:
    """Rotate so the anchor->head vector points along +z. Sign-unambiguous.

    Uses a defined chemical reference (anchor = S, head = chain terminus) so it
    is immune to eigenvector sign flips and near-symmetry degeneracy.
    """
    v = coords[head_idx] - coords[anchor_idx]
    v = v / np.linalg.norm(v)
    target = np.array([0.0, 0.0, 1.0])
    axis = np.cross(v, target)
    n = np.linalg.norm(axis)
    if n < 1e-8:  # already (anti)parallel to z
        rot = np.eye(3) if v[2] > 0 else rotation_matrix("x", 180)
    else:
        axis = axis / n
        angle = math.acos(np.clip(np.dot(v, target), -1.0, 1.0))
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        rot = np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)
    centroid = coords.mean(axis=0)
    return ((coords - centroid) @ rot.T) + centroid
