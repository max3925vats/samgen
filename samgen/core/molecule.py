"""A SAM component: coordinates plus (optional) force-field knowledge.

When an .itp is supplied we parse just enough of it to support anchor detection
and topology assembly: per-atom mass (to find sulfur by element, not by name)
and the bond graph (to tell a terminal/anchor S from an in-chain thioether).
We deliberately do NOT parse the full force field here — topology.py owns that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np

from .gro import GroStructure, read_gro


@dataclass
class Molecule:
    name: str                 # residue name as it appears in the .gro / [molecules]
    struct: GroStructure
    itp_path: Optional[str] = None
    masses: Optional[List[float]] = None          # per-atom, in itp atom order
    bonds: Optional[List[Tuple[int, int]]] = None  # 0-indexed atom pairs

    @property
    def coords(self) -> np.ndarray:
        return self.struct.coords

    @classmethod
    def from_files(cls, name: str, gro: str, itp: Optional[str] = None) -> "Molecule":
        struct = read_gro(gro)
        masses, bonds = (None, None)
        if itp is not None:
            masses, bonds = _parse_itp(itp)
        return cls(name=name, struct=struct, itp_path=itp, masses=masses, bonds=bonds)

    def sulfur_indices(self) -> List[int]:
        """0-indexed atoms whose mass is ~32.06 (sulfur). Name-independent."""
        if self.masses is None:
            raise ValueError(
                f"{self.name}: need an .itp to identify sulfur by mass"
            )
        return [i for i, m in enumerate(self.masses) if abs(m - 32.06) < 0.5]

    def neighbors(self, atom_idx: int) -> List[int]:
        if self.bonds is None:
            raise ValueError(f"{self.name}: need an .itp bond graph")
        out = []
        for a, b in self.bonds:
            if a == atom_idx:
                out.append(b)
            elif b == atom_idx:
                out.append(a)
        return out


def _parse_itp(path: str) -> Tuple[List[float], List[Tuple[int, int]]]:
    """Minimal .itp parse: [atoms] masses + [bonds] graph (first moleculetype)."""
    masses: List[float] = []
    bonds: List[Tuple[int, int]] = []
    section = None
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()  # strip comments
            if not line:
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                continue
            if section == "atoms":
                # nr type resnr resname atomname cgnr charge mass ...
                cols = line.split()
                if len(cols) >= 8:
                    masses.append(float(cols[7]))
            elif section == "bonds":
                cols = line.split()
                if len(cols) >= 2:
                    # itp atom indices are 1-based -> store 0-based
                    bonds.append((int(cols[0]) - 1, int(cols[1]) - 1))
    return masses, bonds
