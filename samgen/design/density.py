"""Density-based / random placement: ligand at a target fraction, seeded.

Deterministic given the seed so surfaces are reproducible. Placement is decided
lazily per site from a per-site hash of (seed, row, col); this avoids needing to
know the lattice extent up front and keeps `label()` pure.
"""

from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass


@dataclass
class Density:
    base: str
    ligand: str
    fraction: float       # target ligand fraction in [0, 1]
    seed: int = 0

    def label(self, row: int, col: int) -> str:
        # Stable uniform value in [0,1) from a hash -> reproducible pattern.
        h = hashlib.sha256(struct.pack("iii", self.seed, row, col)).digest()
        u = int.from_bytes(h[:8], "big") / 2 ** 64
        return self.ligand if u < self.fraction else self.base
