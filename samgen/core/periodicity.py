"""Periodicity checks for a built SAM surface.

Verifies the surface tiles seamlessly under PBC:
  1. the box is a whole multiple of the Au lattice vectors,
  2. the nearest-neighbour strand distance (minimum-image) is >= a threshold
     (catches seams and clashes for any design),
  3. the ligand nearest-neighbour spacing is uniform (no compressed gap at the
     box edge) -- the seam test, stagger-agnostic.

Run on the in-memory structure (full-precision coords) before writing the .gro.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from .gro import GroStructure
from .lattice import Lattice


@dataclass
class PeriodicityReport:
    ok: bool
    issues: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        if self.ok:
            return "periodicity check: OK"
        return "periodicity check FAILED:\n  - " + "\n  - ".join(self.issues)


def _residue_centroids(struct: GroStructure, resname: Optional[str] = None) -> np.ndarray:
    """xy centroid per residue (optionally filtered to `resname`)."""
    groups: dict = {}
    for a in struct.atoms:
        if resname is not None and a.resname != resname:
            continue
        groups.setdefault(a.resid, []).append((a.x, a.y))
    if not groups:
        return np.empty((0, 2))
    return np.array([np.mean(groups[i], axis=0) for i in sorted(groups)])


def _nn_distances(xy: np.ndarray, box: Tuple[float, ...]) -> np.ndarray:
    """Nearest-neighbour distance per point under the minimum-image convention."""
    bx, by = box[0], box[1]
    out = np.empty(len(xy))
    for i in range(len(xy)):
        dx = xy[:, 0] - xy[i, 0]
        dy = xy[:, 1] - xy[i, 1]
        dx -= bx * np.round(dx / bx)
        dy -= by * np.round(dy / by)
        d = np.hypot(dx, dy)
        d[i] = np.inf
        out[i] = d.min()
    return out


def check_box_multiple(box, lat: Lattice, tol: float) -> List[str]:
    issues = []
    for dim, sep, name in ((box[0], lat.colsep, "x"), (box[1], lat.rowsep, "y")):
        n = dim / sep
        if abs(n - round(n)) * sep > tol:
            issues.append(f"box {name}={dim:.4f} nm is not a whole multiple of "
                          f"{sep:.4f} nm (Au lattice does not tile in {name})")
    return issues


def check_min_distance(xy: np.ndarray, box, min_spacing: float) -> List[str]:
    if len(xy) < 2:
        return []
    gmin = float(_nn_distances(xy, box).min())
    if gmin < min_spacing:
        return [f"min strand-strand distance {gmin:.3f} nm < min_spacing "
                f"{min_spacing:.3f} nm (possible seam or clash across PBC)"]
    return []


def check_ligand_uniformity(lig_xy: np.ndarray, box, tol: float) -> List[str]:
    if len(lig_xy) < 2:
        return []
    nn = _nn_distances(lig_xy, box)
    spread = float(nn.max() - nn.min())
    if spread > tol:
        return [f"ligand nearest-neighbour spacing varies by {spread*10:.2f} A "
                f"(> tol {tol*10:.2f} A): ligand grid not uniformly periodic "
                "(seam at the box edge)"]
    return []


def check_surface(struct: GroStructure, lat: Lattice,
                  ligand_resname: Optional[str] = None,
                  min_spacing: float = 0.40, tol: float = 0.002,
                  site_xy: Optional[np.ndarray] = None) -> PeriodicityReport:
    """Run all applicable periodicity checks; return a combined report.

    Parameters
    ----------
    site_xy:
        Optional (N, 2) array of Au placement-site (x, y) coordinates.
        When provided, the min-distance check uses these exact site positions
        rather than per-residue xy centroids.  This avoids false positives on
        mixed surfaces (e.g. coh base + ch3 ligand) where differently-shaped
        strands have centroids that appear closer than the actual Au-site
        separation.  The box-multiple and ligand-uniformity checks are
        unaffected by this parameter.
    """
    issues: List[str] = []
    issues += check_box_multiple(struct.box, lat, tol)
    # Use placement-site coords when available; fall back to residue centroids
    # for surfaces loaded from disk where site info is not stored.
    ref_xy = site_xy if site_xy is not None else _residue_centroids(struct)
    issues += check_min_distance(ref_xy, struct.box, min_spacing)
    if ligand_resname is not None:
        # Uniformity check always uses ligand centroids — correct for a single
        # ligand type (shape is uniform so centroid spacing reflects site spacing).
        issues += check_ligand_uniformity(
            _residue_centroids(struct, ligand_resname), struct.box, tol)
    return PeriodicityReport(ok=not issues, issues=issues)
