"""Minimal GROMACS .top/.itp parser for topology assembly.

Splits a topology file into the pieces we need to merge multiple components into
one simulation-ready topol.top:

    * [defaults]      (kept verbatim; must agree across components)
    * [atomtypes]     (unioned across components, deduped by type name)
    * [moleculetype]  blocks (each molecule's full definition, kept verbatim)

`#include` directives are inlined (resolved relative to the including file) so a
component whose parameters live in a wrapper .top (e.g. a ligand inside a
topol.top that also `#include`s the base .itp) is parsed correctly.
[system]/[molecules] from any source are ignored — we write those ourselves.
"""

from __future__ import annotations

import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

_SECTION_RE = re.compile(r"^\s*\[\s*(\w+)\s*\]")

# Sections that belong to the preceding [moleculetype] block.
_MOL_MEMBER = {
    "atoms", "bonds", "pairs", "pairs_nb", "angles", "dihedrals",
    "exclusions", "constraints", "settles", "virtual_sites2",
    "virtual_sites3", "virtual_sites4", "virtual_sitesn",
    "position_restraints", "distance_restraints", "dihedral_restraints",
    "orientation_restraints",
}


@dataclass
class TopFile:
    defaults: List[str] = field(default_factory=list)         # raw lines incl. header
    atomtypes: "OrderedDict[str, str]" = field(default_factory=OrderedDict)  # name -> raw line
    molecules: "OrderedDict[str, List[str]]" = field(default_factory=OrderedDict)  # name -> raw block

    def atomtype_names_used(self, mol_name: str) -> List[str]:
        """Atom-type names referenced in a molecule's [atoms] (col 2)."""
        used, in_atoms = [], False
        for line in self.molecules.get(mol_name, []):
            m = _SECTION_RE.match(line)
            if m:
                in_atoms = (m.group(1).lower() == "atoms")
                continue
            if in_atoms:
                code = line.split(";", 1)[0].split()
                if len(code) >= 2:
                    used.append(code[1])
        return used


def _inline(path: str, _seen: Optional[set] = None) -> List[str]:
    _seen = _seen or set()
    real = os.path.realpath(path)
    if real in _seen:  # guard against include cycles
        return []
    _seen.add(real)
    base = os.path.dirname(real)
    out: List[str] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            s = line.strip()
            if s.startswith("#include"):
                inc = s.split(None, 1)[1].strip().strip('"<>')
                cand = inc if os.path.isabs(inc) else os.path.join(base, inc)
                if os.path.exists(cand):
                    out.extend(_inline(cand, _seen))
                # silently skip unresolved includes (e.g. forcefield.itp)
            else:
                out.append(line)
    return out


def parse_top(path: str) -> TopFile:
    lines = _inline(path)
    tf = TopFile()
    i, n = 0, len(lines)
    cur_section = None
    cur_mol: Optional[str] = None

    while i < n:
        line = lines[i]
        m = _SECTION_RE.match(line)
        if m:
            cur_section = m.group(1).lower()
            if cur_section == "moleculetype":
                # name appears on the next non-comment, non-blank data line
                name = _next_token(lines, i + 1)
                cur_mol = name
                tf.molecules[name] = [line]
                i += 1
                continue
            if cur_section in ("system", "molecules"):
                cur_mol = None  # stop appending to any molecule
            elif cur_section == "defaults":
                tf.defaults.append(line)
            elif cur_mol and cur_section in _MOL_MEMBER:
                tf.molecules[cur_mol].append(line)
            i += 1
            continue

        # body line
        if cur_section == "defaults":
            tf.defaults.append(line)
        elif cur_section == "atomtypes":
            code = line.split(";", 1)[0].split()
            if code:
                tf.atomtypes.setdefault(code[0], line)
        elif cur_mol and cur_section in (_MOL_MEMBER | {"moleculetype"}):
            tf.molecules[cur_mol].append(line)
        i += 1

    return tf


def _next_token(lines: List[str], start: int) -> str:
    for line in lines[start:]:
        s = line.split(";", 1)[0].strip()
        if s and not s.startswith("["):
            return s.split()[0]
    raise ValueError("moleculetype block has no name line")
