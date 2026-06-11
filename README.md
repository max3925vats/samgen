# samgen

Build self-assembled monolayer (SAM) surfaces as GROMACS `.gro` / `.top` files.

`samgen`'s job is to place strands correctly on a periodic Au(111) grid per a
chosen design, for one-sided and two-sided (shared-sulfur midplane) monolayers.
Orientation/tilt and topology assembly are provided as conveniences that **you
must verify** before running MD.

## What it does

- **Hexagonal Au(111) lattice** — `(√3×√3)R30°` tiling with the staggered-row
  offset, periodic-correct box.
- **SAM orientation** — apply the tilt (`alpha` ≈ 28°, `beta` ≈ 53°; Love et
  al. 2005) with an always-on orientation sanity check and opt-in,
  chemistry-anchored canonicalization.
- **Five placement designs** — density (default), uniform, two-component grid,
  fraction/random (seeded), and multi-ligand (3+).
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

## Placement designs

A design decides which component sits at each Au(111) lattice site. Set it in the
`design:` block of the config.

- **`density`** (default) — ligands at a target areal density (ligands/nm²),
  placed on a uniform sub-lattice of the gold sites. samgen keeps the grid
  perfectly periodic; the realized density is quantized to the Au spacing and is
  reported. If your box and density don't tile uniquely, samgen offers two
  periodic choices (a smaller box with fewer ligands, or a larger box with more)
  and you pick — or set `density_choice: below|above` for batch runs.

  ```yaml
  design: {type: density, base: coh, ligand: ch3, density: 1.0}
  ```

  Example: density 1.0/nm² in a ~10×10 nm box → a 22×24 Au grid, stride 2 →
  **11×12 = 132 ligands** (box 10.978×10.368 nm in rounded-lattice mode).

  For full manual control, `ligand_grid: [nx, ny]` places exactly nx×ny ligands
  (the box is grown so the grid stays periodic; per-axis spacing may differ, so
  anisotropic grids are allowed). It overrides `density` when both are set.

- **`fraction`** — random ligands at a target fraction (0–1), seeded for
  reproducibility.

  ```yaml
  design: {type: fraction, base: coh, ligand: ch3, fraction: 0.25, seed: 42}
  ```

- **`uniform`** — one component on every site.

  ```yaml
  design: {type: uniform, component: ch3}
  ```

- **`grid`** — two components from a 0/1 pattern file.

  ```yaml
  design: {type: grid, pattern: inputs/pattern.dat, mapping: {0: coh, 1: ch3}}
  ```

- **`multilig`** — 3+ components from an integer pattern file.

The config `box: {x, y}` sets the surface area; the lattice fills it and snaps
x,y to whole Au cells for periodicity. The pattern file (for `grid`/`multilig`)
does not set the box — it is sampled per site and padded with the base component
beyond its extent.

### Periodicity check

Every build is checked for periodicity under PBC (on by default, warn-only): the
box is a whole Au-lattice multiple, the ligand grid has no seam (uniform
nearest-neighbour spacing including the wrap), and no two strands fall closer
than `min_spacing` across the periodic boundary. Configure or disable it:

```yaml
periodicity_check: {enabled: true, on_failure: warn, min_spacing: 0.40, tol: 0.002}
```

Set `on_failure: error` to make a non-periodic surface a hard failure instead of
a warning.

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

## Orientation applicability

samgen derives a strand's tilt/twist frame from the anchor sulfur and the first N
backbone carbons of a linear alkyl spacer (default N=9, set per component with
`backbone_carbons`). This matches the canonical alkanethiol-on-Au(111) SAM motif.
Strands that depart from it — non-thiol anchors; branched, short, rigid, or
aromatic spacers; or headgroups attached within the first N carbons — may not be
oriented correctly. In those cases, supply a pre-oriented strand
(`apply_tilt: false`) or verify and correct the generated frame yourself. samgen
does not guarantee orientation for arbitrary strand chemistries.

## Surface types

samgen places strands on a periodic Au(111) grid — that placement is the
guarantee. One-sided and two-sided (shared-sulfur midplane) surfaces are
supported. Orientation and topology assembly are provided as conveniences;
**verify both before running MD** (see the "You finish run-readiness" section and
the "Orientation applicability" section above).

## Two-sided SAMs

The tool builds the **geometry** of a two-sided (shared-S) strand and surface.
Because forming the shared sulfur removes the methyl caps and changes local
charges, you **parameterize the generated strand yourself** (acpype/RESP), then
hand the `.itp`/`.top` back to assemble the integrated topology.
