"""Hexagonal Au(111) lattice for SAM tiling.

Geometry follows Love et al. 2005 (alkanethiolate
(sqrt(3) x sqrt(3))R30 overlayer on Au(111)), lattice constant a = 0.288 nm:

    colsep = sqrt(3) * a   (spacing along a row)
    rowsep = 1.5 * a       (spacing between rows)
    offset = sqrt(3)/2 * a (x-shift applied to alternate rows -> hex packing)

Reference: Love, Estroff, Kriebel, Nuzzo & Whitesides, Chem. Rev. 105 (2005)
1103-1169.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterator, Tuple


@dataclass
class Lattice:
    a: float = 0.288  # Au(111) lattice constant, nm
    # Optional explicit overrides. When None, spacings are computed from `a`.
    # Lattice.rounded() supplies pre-rounded constants for a slightly different
    # (rounded) spacing convention.
    colsep_override: float | None = None
    rowsep_override: float | None = None
    offset_override: float | None = None

    @classmethod
    def rounded(cls) -> "Lattice":
        """Pre-rounded spacing constants for a = 0.288 nm.

        An alternative to the exact sqrt(3)*a spacing: colsep/rowsep/offset are
        rounded to 3-4 decimals. Use when a rounded-constant cell is required.
        """
        return cls(a=0.288, colsep_override=0.499,
                   rowsep_override=0.432, offset_override=0.2494)

    @property
    def colsep(self) -> float:
        return self.colsep_override if self.colsep_override is not None else math.sqrt(3.0) * self.a

    @property
    def rowsep(self) -> float:
        return self.rowsep_override if self.rowsep_override is not None else 1.5 * self.a

    @property
    def offset(self) -> float:
        return self.offset_override if self.offset_override is not None else math.sqrt(3.0) / 2.0 * self.a

    # Tiny epsilon so a strand sitting exactly on the box edge isn't double-
    # counted by the "< box" fill condition.
    _EPS = 1e-9

    def dimensions(self, boxx: float, boxy: float,
                   even_cols: bool = False) -> Tuple[int, int]:
        """Return (ncols, nrows) for a *periodic* tile bounded by boxx x boxy.

        Safeguards keep the final box a valid periodic cell (see docs/DESIGN.md):

        * ncols is fixed by the first (non-offset) row and applied to EVERY row,
          so offset rows aren't one strand short. Offset rows therefore extend
          slightly past boxx by design.
        * nrows is rounded up to an even number, so tiling always ends on a
          complete A/B (non-offset/offset) row pair and stays periodic in y.
        * `even_cols` rounds ncols up to even too. Patterned designs
          (grid/density/multilig) require this so the 2-site stagger lines up
          with the design grid; the uniform design does not.
        """
        ncols = int((boxx - self._EPS) // self.colsep) + 1
        nrows = int((boxy - self._EPS) // self.rowsep) + 1
        if even_cols and ncols % 2 == 1:
            ncols += 1
        if nrows % 2 == 1:          # complete the A/B row pair
            nrows += 1
        return ncols, nrows

    def sites(self, boxx: float, boxy: float,
              even_cols: bool = False) -> Iterator[Tuple[int, int, float, float]]:
        """Yield (row, col, x, y) for every lattice site of the periodic tile.

        Alternate rows are x-shifted by `offset` for hexagonal packing. Every
        row has the same column count (see `dimensions`), so the surface edges
        line up and the reported final box is genuinely periodic.
        """
        ncols, nrows = self.dimensions(boxx, boxy, even_cols)
        for row in range(nrows):
            xstart = self.offset if (row % 2 == 1) else 0.0
            y = row * self.rowsep
            for col in range(ncols):
                yield row, col, xstart + col * self.colsep, y

    def final_box(self, ncols: int, nrows: int, boxz: float) -> Tuple[float, float, float]:
        """Periodic box that matches how many sites were actually placed.

        Slightly different from the requested box: x and y snap to whole
        multiples of colsep/rowsep so inter-strand spacing is preserved.
        """
        return (ncols * self.colsep, nrows * self.rowsep, boxz)
