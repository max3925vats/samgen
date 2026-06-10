# samgen

Build self-assembled monolayer (SAM) surfaces as GROMACS `.gro` / `.top` files.

`samgen` tiles oriented thiol strands onto a hexagonal Au(111) lattice and
assembles a simulation-ready topology, for one-sided and two-sided
(shared-sulfur midplane) monolayers. See [`docs/DESIGN.md`](docs/DESIGN.md) for
the full specification and rationale.

## What it does

- **Hexagonal Au(111) lattice** — `(√3×√3)R30°` tiling with the staggered-row
  offset, periodic-correct box.
- **SAM orientation** — apply the tilt (`alpha` ≈ 28°, `beta` ≈ 53°; Love et
  al. 2005) with an always-on orientation sanity check and opt-in,
  chemistry-anchored canonicalization.
- **Four placement designs** — uniform, two-component grid, density/random
  (seeded), and multi-ligand (3+).
- **Anchor / cap detection** — prompt-first with consent-to-guess.
- **Topology assembly** — residue reordering, merged `[defaults]`/`[atomtypes]`,
  per-component bonded params, and a `[molecules]` block with correct counts,
  optionally gated by `gmx grompp`.
- **Two-sided SAMs** — fuse two arms onto a single shared sulfur via a proper,
  chirality-preserving 180° rotation; emits geometry for you to parameterize.

## Install

```bash
pip install -e .            # needs Python 3.11+, numpy, pyyaml
brew install gromacs        # for boxing/centering and topology validation
```

## Three independent stages

```bash
samgen geometry configs/ch3_onesided.yaml -o sam.gro      # tile (no force field)
samgen topology configs/ch3_grid.yaml --gro sam.gro -o topol.top --validate
samgen build    configs/ch3_grid.yaml -o sam.gro --top topol.top
samgen twosided configs/ch3_grid.yaml --component ch3 -o strand.gro --surface twosided-sam.gro
```

`--validate` runs the `gmx grompp` gate. For `twosided`, the anchor is
prompt-first: it uses the config's `anchor:` if set, otherwise prompts in an
interactive terminal and only auto-detects with your consent. In batch (no TTY)
it requires `anchor:` in the config unless `allow_anchor_autodetect: true`.

Or from Python:

```python
from samgen import generate_geometry, assemble_topology, build
```

## Input files

The examples ship with ready-to-run strands under `configs/inputs/`: a
methyl-terminated (`ch3`) and a hydroxyl-terminated (`coh`) alkanethiol, each as
a single-strand `.gro` (pre-oriented along z) plus a self-contained `.itp`, and
a design `pattern.dat` for the grid. So `samgen build configs/ch3_onesided.yaml`
and `configs/ch3_grid.yaml` run out of the box.

To use **your own** molecules, drop a single-strand `.gro` (canonically oriented
along z, or pre-tilted with `apply_tilt: false`) and a matching `.itp`/`.top`
into `configs/inputs/` and point the config's component paths at them.

## You finish run-readiness

`samgen` builds geometry and a draft topology but deliberately leaves two
simulation-setup steps to you (it prints a reminder when you assemble a
topology):

1. **Position restraints** — restrain the anchor sulfur (and typically the 7th
   chain carbon) so the SAM stays anchored. Add a `[ position_restraints ]`
   block / `#ifdef POSRES` to each component `.itp` and define `POSRES` in your
   `.mdp`.
2. **Center / box the structure** — e.g.
   `gmx editconf -f sam.gro -o centered.gro -c -d <pad>` before
   solvation/equilibration.

## Surface types

Single-component, two-component by design grid, density-based/random, and
multi-ligand (3+). Plus one-sided and two-sided (shared-sulfur midplane) SAMs.

## Two-sided SAMs

The tool builds the **geometry** of a two-sided (shared-S) strand and surface.
Because forming the shared sulfur removes the methyl caps and changes local
charges, you **parameterize the generated strand yourself** (acpype/RESP), then
hand the `.itp`/`.top` back to assemble the integrated topology. See
[`docs/DESIGN.md`](docs/DESIGN.md) §6.
