"""Design strategies: decide which component label occupies each lattice site.

Every strategy implements the same interface so the tiling loop in geometry.py
stays identical across all four surface types:

    Design.label(row, col) -> str   # component key, e.g. "base" or "ligand"

Use `make_design(config)` to construct the right one from a config dict.
"""

from __future__ import annotations

from typing import Dict, Any

from .uniform import Uniform
from .grid import Grid
from .density import Density
from .multilig import MultiLigand


def make_design(config: Dict[str, Any]):
    """Factory: build a Design from the `design:` block of a surface config."""
    dtype = config.get("type", "uniform")
    if dtype == "uniform":
        return Uniform(component=config["component"])
    if dtype == "grid":
        return Grid.from_file(config["pattern"], mapping=config["mapping"])
    if dtype == "density":
        return Density(
            base=config["base"],
            ligand=config["ligand"],
            fraction=config["fraction"],
            seed=config.get("seed", 0),
        )
    if dtype == "multilig":
        return MultiLigand.from_file(config["pattern"], mapping=config["mapping"])
    raise ValueError(f"unknown design type {dtype!r}")


__all__ = ["make_design", "Uniform", "Grid", "Density", "MultiLigand"]
