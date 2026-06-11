"""Command-line interface: `samgen geometry|topology|build CONFIG`.

Mirrors the three independently runnable stages. Anchor handling follows the
prompt-first / consent-to-guess policy: prompt in interactive mode, never
silently auto-detect in batch mode without explicit consent.
"""

from __future__ import annotations

import argparse
import os
import sys

try:
    import yaml
except ImportError:  # keep import-light; yaml only needed for the CLI
    yaml = None

from ._build import build, load_components
from .geometry import generate_geometry, generate_twosided
from .topology import assemble_topology
from .interactive import resolve_anchor_interactive


def _load_config(path: str) -> dict:
    if yaml is None:
        sys.exit("PyYAML is required for the CLI: pip install pyyaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def cmd_geometry(args):
    cfg = _load_config(args.config)
    root = os.path.dirname(os.path.abspath(args.config))
    comps = load_components(cfg, root)
    res = generate_geometry(cfg, comps, out_gro=args.out,
                            manifest_path=args.out + ".manifest.json", root=root)
    print(f"wrote {res.surface_gro} ({res.manifest['natoms']} atoms)")


def cmd_topology(args):
    cfg = _load_config(args.config)
    root = os.path.dirname(os.path.abspath(args.config))
    comps = load_components(cfg, root)
    order = cfg.get("output", {}).get("order") or [m.name for m in comps.values()]
    itp_map = {m.name: m.itp_path for m in comps.values()}
    counts = assemble_topology(args.gro, itp_map=itp_map, order=order,
                               out_top=args.out, out_gro=args.reordered,
                               validate=args.validate)
    print(f"wrote {args.out}: {counts}")


def cmd_build(args):
    cfg = _load_config(args.config)
    root = os.path.dirname(os.path.abspath(args.config))
    gro, top, counts = build(cfg, root, out_gro=args.out, out_top=args.top)
    print(f"wrote {gro} and {top}: {counts}")


def cmd_twosided(args):
    cfg = _load_config(args.config)
    root = os.path.dirname(os.path.abspath(args.config))
    comps = load_components(cfg, root)
    key = args.component or next(iter(comps))
    if key not in comps:
        sys.exit(f"component {key!r} not in config (have: {list(comps)})")
    mol = comps[key]
    specified = cfg["components"][key].get("anchor")
    allow = cfg.get("allow_anchor_autodetect", False)
    res = resolve_anchor_interactive(mol, specified, allow)
    if res.cap_carbon_idx is None:
        sys.exit(f"anchor {mol.struct.atoms[res.anchor_idx].atomname!r} has no "
                 "methyl cap to strip; cannot build a shared-S two-sided strand")
    tsr, geom = generate_twosided(mol, res.anchor_idx, res.cap_carbon_idx, cfg,
                                  out_strand=args.out, out_surface=args.surface)
    print(f"wrote two-sided strand {tsr.strand_gro} ({tsr.natoms} atoms) and "
          f"geometry-only surface {args.surface} ({geom.manifest['natoms']} atoms). "
          "Parameterize the strand, then build the surface from its .itp.")


def main(argv=None):
    p = argparse.ArgumentParser(prog="samgen", description="Build SAM surfaces for GROMACS")
    sub = p.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("geometry", help="tile a surface (no force field)")
    g.add_argument("config")
    g.add_argument("-o", "--out", default="sam.gro")
    g.set_defaults(func=cmd_geometry)

    t = sub.add_parser("topology", help="assemble topol.top for an existing .gro")
    t.add_argument("config")
    t.add_argument("--gro", required=True, help="existing surface .gro")
    t.add_argument("-o", "--out", default="topol.top")
    t.add_argument("--reordered", default="sam-reordered.gro")
    t.add_argument("--validate", action="store_true", help="run gmx grompp gate")
    t.set_defaults(func=cmd_topology)

    b = sub.add_parser("build", help="geometry + topology in one shot")
    b.add_argument("config")
    b.add_argument("-o", "--out", default="sam.gro")
    b.add_argument("--top", default="topol.top")
    b.set_defaults(func=cmd_build)

    w = sub.add_parser("twosided", help="build a two-sided (shared-S) strand + SAM geometry")
    w.add_argument("config")
    w.add_argument("--component", help="which component to fuse (default: first)")
    w.add_argument("-o", "--out", default="twosided-strand.gro",
                   help="output strand .gro (parameterize this)")
    w.add_argument("--surface", default="twosided-sam.gro",
                   help="output full two-sided SAM geometry (geometry-only)")
    w.set_defaults(func=cmd_twosided)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
