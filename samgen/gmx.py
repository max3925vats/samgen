"""Thin wrappers around GROMACS (installed via Homebrew).

We shell out to `gmx` for box/centering (`editconf`) and for topology
validation (`grompp`). Orientation does NOT depend on GROMACS; the
`editconf -princ` cross-check used in tests lives here too.
"""

from __future__ import annotations

import shutil
import subprocess
from typing import Optional


def gmx_binary() -> str:
    """Locate the gmx executable, raising a clear error if missing."""
    for name in ("gmx", "gmx_mpi", "gmx_d"):
        if shutil.which(name):
            return name
    raise RuntimeError(
        "GROMACS not found on PATH. Install with `brew install gromacs`."
    )


def available() -> bool:
    try:
        gmx_binary()
        return True
    except RuntimeError:
        return False


def editconf(in_gro: str, out_gro: str, box: Optional[tuple] = None,
             center: bool = True) -> None:
    """Box and/or center a structure with `gmx editconf`."""
    cmd = [gmx_binary(), "editconf", "-f", in_gro, "-o", out_gro]
    if box is not None:
        cmd += ["-box", *(str(v) for v in box)]
    if center:
        cmd += ["-c"]
    _run(cmd)


def grompp(top: str, gro: str, mdp: str, out_tpr: str = "validate.tpr") -> subprocess.CompletedProcess:
    """Validate a topology+structure. Non-zero return signals a bad topology.

    Used as the acceptance gate for assembled/synthesized topologies: grompp
    hard-fails on missing bonded params, atom-count mismatch, bad exclusions,
    and warns on non-integer charge.
    """
    cmd = [gmx_binary(), "grompp", "-f", mdp, "-c", gro, "-p", top, "-o", out_tpr]
    return _run(cmd, check=False)


def _run(cmd, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=check)
