"""Robust GROMACS .gro reader/writer.

The .gro format is fixed-width by column. We centralise that layout here so the
rest of the package never touches raw columns.

.gro fixed layout (1-indexed columns), positions in nm:
    1-5    residue number
    6-10   residue name
    11-15  atom name
    16-20  atom number
    21-28  x   (%8.3f)
    29-36  y
    37-44  z
    [45-52 53-60 61-68]  optional vx vy vz (%8.4f)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import numpy as np


@dataclass
class GroAtom:
    resid: int
    resname: str
    atomname: str
    atomnum: int
    x: float
    y: float
    z: float
    # Velocities are preserved if present so we never silently drop data.
    vel: Optional[Tuple[float, float, float]] = None


@dataclass
class GroStructure:
    title: str
    atoms: List[GroAtom]
    box: Tuple[float, ...]  # 3 (or 9, for triclinic) box values, nm

    @property
    def coords(self) -> np.ndarray:
        """N x 3 array of positions in nm."""
        return np.array([[a.x, a.y, a.z] for a in self.atoms], dtype=float)

    def set_coords(self, coords: np.ndarray) -> None:
        for atom, (x, y, z) in zip(self.atoms, coords):
            atom.x, atom.y, atom.z = float(x), float(y), float(z)

    @property
    def natoms(self) -> int:
        return len(self.atoms)


def read_gro(path: str) -> GroStructure:
    with open(path) as fh:
        lines = fh.read().splitlines()

    title = lines[0]
    natoms = int(lines[1].strip())
    atoms: List[GroAtom] = []

    for line in lines[2 : 2 + natoms]:
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnum = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        vel = None
        if len(line) >= 68:  # velocities present
            vel = (float(line[44:52]), float(line[52:60]), float(line[60:68]))
        atoms.append(GroAtom(resid, resname, atomname, atomnum, x, y, z, vel))

    box = tuple(float(v) for v in lines[2 + natoms].split())
    return GroStructure(title=title, atoms=atoms, box=box)


def write_gro(struct: GroStructure, path: str) -> None:
    lines = [struct.title, f"{struct.natoms:5d}"]
    for a in struct.atoms:
        # GROMACS truncates resid/atomnum to 5 digits (mod 100000); match it so
        # large surfaces stay valid .gro files.
        resid = a.resid % 100000
        atomnum = a.atomnum % 100000
        line = (
            f"{resid:5d}{a.resname:<5.5s}{a.atomname:>5.5s}{atomnum:5d}"
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
        )
        if a.vel is not None:
            line += f"{a.vel[0]:8.4f}{a.vel[1]:8.4f}{a.vel[2]:8.4f}"
        lines.append(line)
    lines.append(" ".join(f"{v:.5f}" for v in struct.box))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
