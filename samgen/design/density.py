"""Areal-density placement on a uniform Au(111) sub-lattice.

Given a target areal ligand density (ligands/nm^2), place ligands on a uniform
sub-lattice of the Au sites chosen so the result is ALWAYS perfectly periodic
(no compressed gap at the box boundary). The realized density is quantized to
rho_full / k^2 for an integer stride k; when the requested box does not tile the
stride, two periodic options bracketing the box are offered (see the design class
and interactive resolution).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple


def choose_stride(site_density: float, target: float) -> int:
    """Isotropic stride k whose realized density (site_density/k^2) is nearest target."""
    if target <= 0:
        raise ValueError("density must be > 0 ligands/nm^2")
    return max(1, round(math.sqrt(site_density / target)))


def _axis_multiples(n: int, k: int) -> Tuple[int, int]:
    """Lower and upper multiples of k bracketing n (equal when k divides n)."""
    if n % k == 0:
        return (n, n)
    lo = (n // k) * k
    return (lo, lo + k)


@dataclass(frozen=True)
class GridOption:
    ncols: int          # Au columns in the (possibly grown) box
    nrows: int          # Au rows
    kx: int             # column stride (ligand every kx-th column)
    ky: int             # row stride
    box_x: float
    box_y: float

    @property
    def nx(self) -> int:
        return self.ncols // self.kx

    @property
    def ny(self) -> int:
        return self.nrows // self.ky

    @property
    def count(self) -> int:
        return self.nx * self.ny

    @property
    def density(self) -> float:
        return self.count / (self.box_x * self.box_y)


def grid_options(ncols: int, nrows: int, k: int,
                 colsep: float, rowsep: float) -> List[GridOption]:
    """Perfectly-periodic isotropic ligand grids near (ncols, nrows) at stride k.

    Returns 1 option when k divides both axes (unique), else up to 4 (the
    box bracketed to lower/upper multiples of k on each indivisible axis).
    Both strides equal k (isotropic).
    """
    xs = sorted(set(_axis_multiples(ncols, k)))
    ys = sorted(set(_axis_multiples(nrows, k)))
    # Filter out degenerate options where an axis multiple rounded down to zero
    # (happens when k > ncols or k > nrows at extreme low density).
    return [
        GridOption(cx, cy, k, k, cx * colsep, cy * rowsep)
        for cx in xs for cy in ys
        if cx > 0 and cy > 0
    ]


def explicit_grid(ncols: int, nrows: int, nx: int, ny: int,
                  colsep: float, rowsep: float) -> GridOption:
    """An explicit nx x ny ligand grid (ligand_grid override).

    Per-axis strides are chosen independently (anisotropic allowed) and the box
    is grown/snapped so the grid is perfectly periodic: ncols = nx*kx, etc.
    """
    if nx < 1 or ny < 1:
        raise ValueError("ligand_grid entries must be >= 1")
    kx = max(1, round(ncols / nx))
    ky = max(1, round(nrows / ny))
    return GridOption(nx * kx, ny * ky, kx, ky, nx * kx * colsep, ny * ky * rowsep)


def select_grid(options: List[GridOption], choice: Optional[str]) -> Optional[GridOption]:
    """Pick an option without prompting, or return None if a choice is needed.

    Unique -> that option. choice 'below'/'above' -> min/max ligand count.
    """
    if len(options) == 1:
        return options[0]
    if choice == "below":
        return min(options, key=lambda o: o.count)
    if choice == "above":
        return max(options, key=lambda o: o.count)
    return None


class Density:
    """Areal-density design. Needs the full grid, so it is resolved (per-axis
    strides kx, ky and final ncols/nrows chosen) before tiling; see
    interactive.resolve_density_interactive.
    """

    def __init__(self, base: str, ligand: str, density: float,
                 choice: Optional[str] = None,
                 ligand_grid: Optional[Tuple[int, int]] = None) -> None:
        self.base = base
        self.ligand = ligand
        self.density = density          # target ligands/nm^2
        self.choice = choice            # 'below' | 'above' for batch mode
        self.ligand_grid = ligand_grid  # explicit (nx, ny) override, or None
        self._kx = 0                    # set by configure()
        self._ky = 0

    def configure(self, kx: int, ky: int) -> None:
        self._kx, self._ky = kx, ky

    def label(self, row: int, col: int) -> str:
        if self._kx <= 0 or self._ky <= 0:
            raise ValueError("Density design not resolved; call configure(kx, ky) first")
        return (self.ligand if (row % self._ky == 0 and col % self._kx == 0)
                else self.base)
