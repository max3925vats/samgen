"""samgen — build self-assembled monolayer (SAM) surfaces for GROMACS.

Three independently runnable stages (see docs/DESIGN.md):
    generate_geometry()  -> tiled surface .gro (no force field needed)
    assemble_topology()  -> integrated topol.top + reordered .gro
    build()              -> geometry then topology in one shot
"""

from .geometry import generate_geometry, build_twosided_strand, generate_twosided
from .topology import assemble_topology
from ._build import build

__all__ = [
    "generate_geometry", "build_twosided_strand", "generate_twosided",
    "assemble_topology", "build",
]
__version__ = "0.0.1"
