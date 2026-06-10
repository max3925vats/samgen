"""Stage 1: geometry generation (no force field required).

Builds a tiled SAM surface .gro from oriented strands. Also writes a sidecar
manifest so the topology stage can run independently without re-deriving counts
or residue order.

Covers one-sided tiling and two-sided strand fusion + full two-sided SAM
geometry (proper 180-deg rotation about a horizontal axis through the shared S;
see docs/DESIGN.md sec. 6).
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, Optional
import numpy as np

from .core.gro import GroStructure, GroAtom, write_gro
from .core.lattice import Lattice
from .core.molecule import Molecule
from .core import orient
from .design import make_design
from .design.density import Density
from .interactive import resolve_density_interactive


@dataclass
class GeometryResult:
    surface_gro: str
    manifest: dict
    twosided_strand_gro: Optional[str] = None
    twosided_surface_gro: Optional[str] = None


def _make_lattice(latcfg: dict) -> Lattice:
    """Build a Lattice from config. `rounded: true` uses pre-rounded spacings."""
    if latcfg.get("rounded"):
        return Lattice.rounded()
    return Lattice(a=latcfg.get("a", 0.288),
                   colsep_override=latcfg.get("colsep"),
                   rowsep_override=latcfg.get("rowsep"),
                   offset_override=latcfg.get("offset"))


def generate_geometry(config: dict, components: Dict[str, Molecule],
                      out_gro: str, manifest_path: Optional[str] = None,
                      root: str = ".",
                      input_fn=input, is_tty=None) -> GeometryResult:
    """Tile `components` onto the lattice per `config`, write `out_gro`.

    `components` maps the keys used by the design (e.g. "base", "ligand") to
    already-loaded Molecules. By default it applies the SAM tilt to each
    (canonical, along-z) strand. Set config["apply_tilt"] = False when the input
    strands are ALREADY oriented and tilted so tiling just places them verbatim.

    `root` resolves a relative design `pattern:` path (grid/multilig) the same
    way component .gro/.itp paths are resolved.
    """
    lat = _make_lattice(config["lattice"])
    alpha = config["lattice"].get("tilt_alpha", 0)
    beta = config["lattice"].get("tilt_beta", 0)
    boxx, boxy, boxz = config["box"]["x"], config["box"]["y"], config["box"]["z"]
    apply_tilt = config.get("apply_tilt", True)

    # Prepare per-component coordinates. When applying the tilt we first run the
    # loud orientation sanity check; when trusting pre-tilted input we skip it.
    tilted: Dict[str, np.ndarray] = {}
    for key, mol in components.items():
        if apply_tilt:
            orient.check_oriented(mol.coords)  # raises if mis-oriented
            tilted[key] = orient.apply_tilt(mol.coords, alpha, beta)
        else:
            tilted[key] = mol.coords

    # Resolve a relative pattern path against the config directory before the
    # design module opens it.
    design_cfg = dict(config["design"])
    if design_cfg.get("pattern") and not os.path.isabs(design_cfg["pattern"]):
        design_cfg["pattern"] = os.path.join(root, design_cfg["pattern"])
    design = make_design(design_cfg)

    # Patterned designs force an even column count; uniform does not.
    even_cols = config["design"].get("type", "density") != "uniform"
    ncols, nrows = lat.dimensions(boxx, boxy, even_cols)

    # The density design resolves a perfectly-periodic stride and may grow the
    # box; this can override (ncols, nrows).
    if isinstance(design, Density):
        opt = resolve_density_interactive(design, lat, ncols, nrows,
                                          input_fn=input_fn, is_tty=is_tty)
        design.configure(opt.kx, opt.ky)
        ncols, nrows = opt.ncols, opt.nrows

    atoms = []
    resid = 0
    counts: Dict[str, int] = {key: 0 for key in components}

    for row, col, x, y in lat.sites_for(ncols, nrows):
        key = design.label(row, col)
        mol, coords = components[key], tilted[key]
        resid += 1
        counts[key] += 1
        shift = np.array([x, y, 0.0])
        for i, atom in enumerate(mol.struct.atoms):
            px, py, pz = coords[i] + shift
            atoms.append(GroAtom(resid, atom.resname, atom.atomname,
                                 len(atoms) + 1, px, py, pz))

    box = lat.final_box(ncols, nrows, boxz)
    struct = GroStructure(title="SAM surface (samgen geometry)", atoms=atoms, box=box)
    write_gro(struct, out_gro)

    manifest = {
        "natoms": struct.natoms,
        "requested_box": [boxx, boxy, boxz],
        "final_box": [round(v, 5) for v in box],  # snapped to lattice multiples
        "grid": {"ncols": ncols, "nrows": nrows},
        "counts": counts,                       # residues per component key
        "resnames": {k: components[k].name for k in components},
        "itp": {components[k].name: components[k].itp_path for k in components},
        "order": config.get("output", {}).get("order"),
        "leaflets": config.get("leaflets", 1),
    }
    if manifest_path:
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)

    return GeometryResult(surface_gro=out_gro, manifest=manifest)


@dataclass
class TwoSidedResult:
    strand: Molecule        # the fused two-sided strand (in-memory)
    strand_gro: str         # path written for the user to parameterize
    natoms: int
    n_armA: int             # atoms in arm A (includes the shared S)
    n_armB: int             # atoms in arm B (no S; mirror of arm A's chain)
    s_index: int            # 0-based index of the shared sulfur in the output


def build_twosided_strand(strand: Molecule, anchor_idx: int,
                          cap_carbon_idx: int, out_gro: str,
                          resname: Optional[str] = None,
                          armB_suffix: str = "b") -> TwoSidedResult:
    """Fuse a one-sided strand into a two-sided (shared-S) strand. GEOMETRY ONLY.

    Approach (docs/DESIGN.md sec. 6):
      1. Canonicalize: anchor->head along +z, shared S translated to the origin.
      2. Strip the methyl cap on the anchor S (cap carbon + its hydrogens).
      3. Arm B = a copy of arm A's non-sulfur atoms, transformed by a PROPER
         180-deg rotation about the x-axis through S. Proper (det +1) so it
         preserves chirality and the genuine beta twist (a reflection would
         build the enantiomer); it flips the chain to point along -z so the two
         arms are antiparallel about the shared S.
      4. Fuse: arm A (with S) + arm B (no S). The result is antiparallel along
         z so the tiling tilt later gives BOTH arms the correct SAM tilt.

    The output .gro is for the user to parameterize externally (the caps are
    gone and charges must be refit); samgen never synthesizes parameters.
    """
    if strand.bonds is None:
        raise ValueError("two-sided fusion needs the strand's bond graph (.itp)")

    # 1. canonicalize: chain along +z, then put the shared S at the origin.
    head_idx = int(np.argmax(np.linalg.norm(strand.coords - strand.coords[anchor_idx], axis=1)))
    coords = orient.canonicalize(strand.coords, anchor_idx, head_idx)
    coords = coords - coords[anchor_idx]

    # 2. cap atoms to strip: the cap carbon + the hydrogens bonded to it.
    cap = {cap_carbon_idx}
    for nb in strand.neighbors(cap_carbon_idx):
        if strand.masses is not None and strand.masses[nb] < 11.0:
            cap.add(nb)

    armA_idx = [i for i in range(len(strand.struct.atoms)) if i not in cap]
    armB_src = [i for i in armA_idx if i != anchor_idx]   # arm B has no sulfur

    # 3. arm B coords via proper 180-deg rotation about x through S (origin).
    Rx = orient.rotation_matrix("x", 180.0)
    coordsB = coords @ Rx.T

    # 4. assemble output atoms: arm A (with S) then arm B.
    name = resname or strand.name
    atoms = []
    s_index = None
    for i in armA_idx:
        a = strand.struct.atoms[i]
        if i == anchor_idx:
            s_index = len(atoms)
        atoms.append(GroAtom(1, name, a.atomname, len(atoms) + 1,
                             *coords[i]))
    for i in armB_src:
        a = strand.struct.atoms[i]
        nm = (a.atomname + armB_suffix)[:5] if len(a.atomname) < 5 else a.atomname
        atoms.append(GroAtom(1, name, nm, len(atoms) + 1, *coordsB[i]))

    box = _padded_box(np.array([[a.x, a.y, a.z] for a in atoms]))
    struct = GroStructure(
        title="two-sided SAM strand (samgen GEOMETRY ONLY - parameterize before use)",
        atoms=atoms, box=box)
    write_gro(struct, out_gro)

    fused = Molecule(name=name, struct=struct, itp_path=None)
    return TwoSidedResult(strand=fused, strand_gro=out_gro, natoms=len(atoms),
                          n_armA=len(armA_idx), n_armB=len(armB_src), s_index=s_index)


def _padded_box(coords: np.ndarray, pad: float = 1.0):
    span = coords.max(axis=0) - coords.min(axis=0)
    return tuple(float(s + 2 * pad) for s in span)


def generate_twosided(strand: Molecule, anchor_idx: int, cap_carbon_idx: int,
                      config: dict, out_strand: str, out_surface: str):
    """Build the two-sided strand AND the full two-sided SAM geometry.

    Both outputs are GEOMETRY ONLY. Parameterize `out_strand`, then build the
    simulation-ready surface from the returned .itp via the normal build path.
    """
    res = build_twosided_strand(strand, anchor_idx, cap_carbon_idx, out_strand)
    cfg = {**config, "design": {"type": "uniform", "component": "twosided"}}
    geom = generate_geometry(cfg, {"twosided": res.strand}, out_gro=out_surface)
    return res, geom
