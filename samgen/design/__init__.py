"""Design strategies: decide which component label occupies each lattice site.

Per-site designs implement `label(row, col) -> str`. The `density` design also
needs the full grid and is resolved before tiling (see geometry/interactive).

Use `make_design(config)` to construct the right one from a config dict.
"""

from __future__ import annotations

from typing import Dict, Any

from .uniform import Uniform
from .grid import Grid
from .fraction import Fraction
from .density import Density
from .multilig import MultiLigand


def make_design(config: Dict[str, Any]):
    """Factory: build a Design from the `design:` block of a surface config.

    Default type is `density` (areal ligand density on a uniform Au sub-lattice).
    """
    dtype = config.get("type", "density")
    if dtype == "uniform":
        return Uniform(component=config["component"])
    if dtype == "grid":
        return Grid.from_file(config["pattern"], mapping=config["mapping"])
    if dtype == "fraction":
        return Fraction(base=config["base"], ligand=config["ligand"],
                        fraction=config["fraction"], seed=config.get("seed", 0))
    if dtype == "density":
        lg = config.get("ligand_grid")
        return Density(base=config["base"], ligand=config["ligand"],
                       density=config.get("density"),
                       choice=config.get("density_choice"),
                       ligand_grid=tuple(lg) if lg else None)
    if dtype == "multilig":
        return MultiLigand.from_file(config["pattern"], mapping=config["mapping"])
    raise ValueError(f"unknown design type {dtype!r}")


__all__ = ["make_design", "Uniform", "Grid", "Fraction", "Density", "MultiLigand"]
