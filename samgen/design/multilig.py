"""Multi-ligand (3+) placement from an integer pattern file.

Identical mechanics to Grid but the mapping carries 3+ component keys, e.g.
{0: "base", 1: "ligA", 2: "ligB", 3: "ligC"}.
"""

from __future__ import annotations

from .grid import Grid


class MultiLigand(Grid):
    """Same pattern semantics as Grid; separate type for clarity/validation."""
    pass
