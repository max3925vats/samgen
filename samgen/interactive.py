"""Interactive anchor resolution: prompt first, guess only on consent.

Policy:
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
from .core.lattice import Lattice
from .core.molecule import Molecule
from .design import density as _density
from .design.density import Density, GridOption


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


def _format_density_options(options, target: float) -> str:
    lines = [f"  target density {target:.3f}/nm^2 - choose a perfectly periodic grid:"]
    for o in sorted(options, key=lambda x: x.count):
        tag = "below" if o.count == min(x.count for x in options) else "above"
        lines.append(f"    [{tag}] {o.nx} x {o.ny} = {o.count} ligands, "
                     f"box {o.box_x:.3f} x {o.box_y:.3f} nm, "
                     f"density {o.density:.3f}/nm^2")
    return "\n".join(lines)


def resolve_density_interactive(
    design: Density,
    lat: Lattice,
    ncols0: int,
    nrows0: int,
    input_fn: Callable[[str], str] = input,
    is_tty: Optional[bool] = None,
) -> GridOption:
    """Resolve a Density design to one perfectly-periodic GridOption.

    Explicit ligand_grid -> built directly (no prompt). Unique stride grid ->
    returned directly. Otherwise: batch mode requires design.choice
    ('below'/'above') or fails with both options printed; interactive mode prompts.
    """
    if design.ligand_grid is not None:
        nx, ny = design.ligand_grid
        return _density.explicit_grid(ncols0, nrows0, nx, ny, lat.colsep, lat.rowsep)

    if design.density is None:
        raise ValueError(
            "density design requires a 'density:' (ligands/nm^2) or a "
            "'ligand_grid:' [nx, ny]"
        )

    k = _density.choose_stride(lat.site_density(), design.density)
    options = _density.grid_options(ncols0, nrows0, k, lat.colsep, lat.rowsep)
    chosen = _density.select_grid(options, design.choice)
    if chosen is not None:
        return chosen

    if is_tty is None:
        is_tty = sys.stdin.isatty()
    if not is_tty:
        raise ValueError(
            f"density {design.density:.3f}/nm^2 does not tile this box uniquely; "
            f"set design.density_choice to 'below' or 'above'.\n"
            + _format_density_options(options, design.density)
        )

    print(_format_density_options(options, design.density))
    ans = input_fn("Pick [below/above]: ").strip().lower()
    if ans not in ("below", "above"):
        ans = "above"
    return _density.select_grid(options, ans)
