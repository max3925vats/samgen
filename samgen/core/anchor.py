"""Anchor (gold-facing S) and methyl-cap detection.

Policy: prompt the user first; auto-detect ONLY with explicit consent, and
always report what was picked and why. Detection never silently guesses.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional
import numpy as np

from .molecule import Molecule

# Carbon mass (GAFF C2 united-atom CH2 ~14.03, c3 ~12.01). A methyl cap carbon
# in these models is a carbon bonded to the anchor S and to nothing else heavy.
_H_MASS = (1.0, 1.1)


@dataclass
class AnchorResult:
    anchor_idx: int           # 0-indexed S atom
    cap_carbon_idx: Optional[int]   # methyl cap carbon bonded to the anchor S
    reason: str               # human-readable explanation (for the consent prompt)


def resolve_anchor(
    mol: Molecule,
    specified: Optional[str | int],
    allow_autodetect: bool,
) -> AnchorResult:
    """Resolve the anchor atom.

    `specified` may be an atom name (e.g. "S41") or a 1-based index. If None and
    autodetect is not permitted, raise so batch runs fail clearly.
    """
    if specified is not None:
        idx = _resolve_specified(mol, specified)
        return AnchorResult(idx, _find_cap(mol, idx), f"specified ({specified})")

    if not allow_autodetect:
        raise ValueError(
            f"{mol.name}: no anchor specified. Set it in the config, or pass "
            "allow_anchor_autodetect: true to consent to auto-detection."
        )
    return autodetect_anchor(mol)


def autodetect_anchor(mol: Molecule) -> AnchorResult:
    """Heuristic anchor detection. Only call after the user has consented.

    Anchor = the gold-facing sulfur, detected as the S bonded to a terminal
    methyl cap. Falls back to the lowest-z sulfur. Raises on ambiguity.
    """
    sulfurs = mol.sulfur_indices()
    if not sulfurs:
        raise ValueError(f"{mol.name}: no sulfur found; specify the anchor manually")

    capped = [(s, c) for s in sulfurs if (c := _find_cap(mol, s)) is not None]
    if len(capped) == 1:
        s, c = capped[0]
        name = mol.struct.atoms[s].atomname
        return AnchorResult(s, c, f"S {name!r} bonded to methyl cap "
                                  f"{mol.struct.atoms[c].atomname!r}")
    if len(capped) > 1:
        names = ", ".join(mol.struct.atoms[s].atomname for s, _ in capped)
        raise ValueError(
            f"{mol.name}: ambiguous anchor — multiple capped sulfurs ({names}). "
            "Specify the anchor manually."
        )

    # Fallback: lowest-z sulfur (assumes a pre-oriented, S-down strand).
    coords = mol.coords
    s = min(sulfurs, key=lambda i: coords[i][2])
    name = mol.struct.atoms[s].atomname
    return AnchorResult(s, None, f"lowest-z sulfur {name!r} (no methyl cap found)")


def _resolve_specified(mol: Molecule, specified: str | int) -> int:
    if isinstance(specified, int):
        return specified - 1  # 1-based -> 0-based
    s = str(specified)
    if s.isdigit():
        return int(s) - 1
    for i, atom in enumerate(mol.struct.atoms):
        if atom.atomname == s:
            return i
    raise ValueError(f"{mol.name}: anchor atom {specified!r} not found")


def backbone_head(mol: Molecule, anchor_idx: int, n_carbons: int = 9) -> int:
    """Index of the Nth alkyl backbone carbon from the anchor S.

    Walks the linear carbon chain from the anchor (skipping the methyl cap).
    Used as the orientation 'head' so the tilt/twist axis follows the alkyl
    spacer, not a divergent ligand headgroup. Warns and stops early if the
    chain branches or ends before N carbons.
    """
    if mol.bonds is None or mol.masses is None:
        raise ValueError(f"{mol.name}: backbone detection needs an .itp bond graph")

    cap = _find_cap(mol, anchor_idx)
    # approximation: O/N bonded in-chain would also pass this carbon test
    # (fine for thiol alkyl SAMs).
    starts = [nb for nb in mol.neighbors(anchor_idx)
              if mol.masses[nb] >= 11.0 and nb != cap]
    if not starts:
        raise ValueError(f"{mol.name}: no alkyl chain carbon bonded to the anchor")

    prev, cur, count = anchor_idx, starts[0], 1
    while count < n_carbons:
        nxt = [nb for nb in mol.neighbors(cur)
               if mol.masses[nb] >= 11.0 and nb != prev]
        if len(nxt) == 0:
            warnings.warn(f"{mol.name}: alkyl chain ends after {count} carbons "
                          f"(< {n_carbons}); orienting on the shorter segment")
            break
        if len(nxt) > 1:
            warnings.warn(f"{mol.name}: alkyl chain branches at carbon {count}; "
                          f"orienting on the segment up to the branch")
            break
        prev, cur = cur, nxt[0]
        count += 1
    return cur


def _find_cap(mol: Molecule, sulfur_idx: int) -> Optional[int]:
    """A methyl cap carbon: bonded to the S and to 3 H and no other heavy atom."""
    if mol.bonds is None or mol.masses is None:
        return None
    for nb in mol.neighbors(sulfur_idx):
        if mol.masses[nb] < 11.0:  # not a carbon
            continue
        heavy = [x for x in mol.neighbors(nb) if mol.masses[x] >= 11.0]
        hyd = [x for x in mol.neighbors(nb) if mol.masses[x] < 11.0]
        # carbon's only heavy neighbour is the sulfur, plus ~3 hydrogens
        if heavy == [sulfur_idx] and len(hyd) >= 2:
            return nb
    return None
