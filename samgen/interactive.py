"""Interactive anchor resolution: prompt first, guess only on consent.

Policy (docs/DESIGN.md sec. 5):
  * If the anchor is specified (config or argument), use it without prompting.
  * In an interactive terminal, prompt the user to type the anchor; if they
    leave it blank, ask whether to auto-detect, then SHOW the guess and its
    reason and require confirmation before using it.
  * Non-interactive (batch) runs never prompt: they defer to the config, and
    only auto-detect if `allow_anchor_autodetect` consent was given.

`input_fn` and `is_tty` are injectable so the flow is unit-testable without a
real terminal.
"""

from __future__ import annotations

import sys
from typing import Callable, Optional

from .core.anchor import AnchorResult, resolve_anchor, autodetect_anchor
from .core.molecule import Molecule


def resolve_anchor_interactive(
    mol: Molecule,
    specified: Optional[str | int],
    allow_autodetect: bool,
    input_fn: Callable[[str], str] = input,
    is_tty: Optional[bool] = None,
) -> AnchorResult:
    if specified is not None:
        return resolve_anchor(mol, specified, allow_autodetect)

    if is_tty is None:
        is_tty = sys.stdin.isatty()
    if not is_tty:
        # batch: resolve_anchor raises unless allow_autodetect consent was set
        return resolve_anchor(mol, None, allow_autodetect)

    typed = input_fn(
        f"Anchor atom for '{mol.name}' (name or 1-based index, "
        "blank to auto-detect): "
    ).strip()
    if typed:
        return resolve_anchor(mol, typed, allow_autodetect)

    consent = input_fn(
        "No anchor given. Attempt auto-detection of the anchor sulfur? [y/N]: "
    ).strip().lower()
    if consent not in ("y", "yes"):
        raise ValueError("anchor not provided and auto-detection declined")

    res = autodetect_anchor(mol)
    name = mol.struct.atoms[res.anchor_idx].atomname
    confirm = input_fn(
        f"Detected anchor {name!r} ({res.reason}). Use it? [y/N]: "
    ).strip().lower()
    if confirm not in ("y", "yes"):
        raise ValueError("auto-detected anchor rejected; specify it explicitly")
    return res
