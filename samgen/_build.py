"""Stage 3: `build` = geometry then topology in one shot.

Convenience path for one-sided surfaces (and two-sided where the parameterized
strand already exists). Loads components from the config, runs geometry, then
assembles the topology against the generated surface.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

from .core.molecule import Molecule
from .geometry import generate_geometry
from .topology import assemble_topology


def load_components(config: dict, root: str = ".") -> Dict[str, Molecule]:
    """Instantiate Molecules for every key in config['components']."""
    comps: Dict[str, Molecule] = {}
    for key, spec in config["components"].items():
        gro = os.path.join(root, spec["gro"])
        itp = os.path.join(root, spec["itp"]) if spec.get("itp") else None
        comps[key] = Molecule.from_files(name=spec["resname"], gro=gro, itp=itp)
    return comps


def build(config: dict, root: str = ".", out_gro: str = "sam.gro",
          out_top: str = "topol.top", out_reordered: Optional[str] = "sam-reordered.gro",
          input_fn=input, is_tty=None):
    """Full pipeline. Returns (surface_gro, top, counts)."""
    components = load_components(config, root)
    out = config.get("output", {})
    geom = generate_geometry(config, components, out_gro=out_gro,
                             manifest_path=out_gro + ".manifest.json", root=root,
                             input_fn=input_fn, is_tty=is_tty)

    order = out.get("order") or [m.name for m in components.values()]
    itp_map = {m.name: m.itp_path for m in components.values()}
    counts = assemble_topology(out_gro, itp_map=itp_map, order=order,
                               out_top=out_top, out_gro=out_reordered)
    return geom.surface_gro, out_top, counts
