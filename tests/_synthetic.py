"""Synthetic SAM strands for tests — exercise the math, not a real molecule."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from samgen.core.gro import GroStructure, GroAtom
from samgen.core.molecule import Molecule

_S_MASS, _C_MASS, _H_MASS = 32.06, 14.027, 1.008


class _SynthMolecule(Molecule):
    anchor_index: int
    _backbone: List[int]

    def backbone_index(self, n: int) -> int:
        """0-based atom index of the n-th backbone carbon (1-based n)."""
        return self._backbone[min(n, len(self._backbone)) - 1]

    def unoriented_copy(self) -> "Molecule":
        """A copy with the chain rotated to lie along x (needs canonicalization)."""
        c = self.coords.copy()
        # swap z<->x so the along-z chain now lies along x
        c = c[:, [2, 1, 0]]
        atoms = [GroAtom(a.resid, a.resname, a.atomname, a.atomnum,
                         float(x), float(y), float(z), a.vel)
                 for a, (x, y, z) in zip(self.struct.atoms, c)]
        m = copy.copy(self)
        m.struct = GroStructure(self.struct.title, atoms, self.struct.box)
        return m


def linear_strand(name: str = "LIG", n_chain: int = 11, with_cap: bool = True,
                  headgroup_offset: Optional[Tuple[float, float, float]] = None
                  ) -> _SynthMolecule:
    """Build an along-z thiol strand: S at z=0, n_chain carbons going +z,
    an optional methyl cap below S, an optional sideways headgroup on the top
    carbon. Carries masses + bonds so anchor/backbone logic works.
    """
    atoms: List[GroAtom] = []
    masses: List[float] = []
    bonds: List[Tuple[int, int]] = []
    bond_len = 0.15

    def add(atomname, x, y, z, mass):
        atoms.append(GroAtom(1, name, atomname, len(atoms) + 1, x, y, z))
        masses.append(mass)
        return len(atoms) - 1

    anchor = add("S1", 0.0, 0.0, 0.0, _S_MASS)
    backbone: List[int] = []
    prev = anchor
    for i in range(n_chain):
        ci = add(f"C{i+1}", 0.0, 0.0, (i + 1) * bond_len, _C_MASS)
        bonds.append((prev, ci))
        backbone.append(ci)
        prev = ci
    if headgroup_offset is not None:
        ox, oy, oz = headgroup_offset
        top = backbone[-1]
        hg = add("O1", atoms[top].x + ox, atoms[top].y + oy, atoms[top].z + oz, 16.0)
        bonds.append((top, hg))
    if with_cap:
        cap = add("Cc", 0.0, 0.0, -bond_len, _C_MASS)
        bonds.append((anchor, cap))
        for h in range(3):
            hi = add(f"Hc{h+1}", 0.05 * (h - 1), 0.05, -bond_len - 0.05, _H_MASS)
            bonds.append((cap, hi))

    box = (5.0, 5.0, 5.0)
    struct = GroStructure(f"synthetic {name}", atoms, box)
    mol = _SynthMolecule(name=name, struct=struct, itp_path=None,
                         masses=masses, bonds=bonds)
    mol.anchor_index = anchor
    mol._backbone = backbone
    return mol
