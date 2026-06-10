"""Unit tests for the working core modules."""

import math
import numpy as np

from samgen.core.gro import GroStructure, GroAtom, write_gro, read_gro
from samgen.core.lattice import Lattice
from samgen.core import orient


def test_gro_roundtrip(tmp_path):
    atoms = [
        GroAtom(1, "COH", "C1", 1, 0.123, 0.456, 0.789),
        GroAtom(1, "COH", "S13", 2, -0.239, 0.107, -0.331),
    ]
    s = GroStructure("t", atoms, (4.0, 4.0, 4.0))
    p = tmp_path / "x.gro"
    write_gro(s, str(p))
    back = read_gro(str(p))
    assert back.natoms == 2
    assert back.atoms[1].atomname == "S13"
    assert math.isclose(back.atoms[0].x, 0.123, abs_tol=1e-3)


def test_lattice_geometry():
    lat = Lattice(a=0.288)
    assert math.isclose(lat.colsep, math.sqrt(3) * 0.288)
    assert math.isclose(lat.rowsep, 1.5 * 0.288)
    sites = list(lat.sites(boxx=1.0, boxy=1.0))
    assert len(sites) > 0
    # alternate rows are x-offset (hex packing)
    row0 = [s for s in sites if s[0] == 0]
    row1 = [s for s in sites if s[0] == 1]
    assert math.isclose(row0[0][2], 0.0)
    assert math.isclose(row1[0][2], lat.offset)


def test_lattice_periodic_equal_rows():
    # Every row must have the same column count, and total rows must be even
    # (complete A/B pairs) so the final box is a valid periodic cell.
    lat = Lattice(a=0.288)
    sites = list(lat.sites(boxx=5.0, boxy=5.0))
    from collections import Counter
    per_row = Counter(r for r, c, x, y in sites)
    assert len(set(per_row.values())) == 1, "rows have unequal strand counts"
    ncols, nrows = lat.dimensions(5.0, 5.0)
    assert nrows % 2 == 0
    assert ncols == 11 and nrows == 12  # expected hex tiling for a 5x5 nm box
    assert len(sites) == ncols * nrows


def test_final_box_snaps_and_differs_from_request():
    lat = Lattice(a=0.288)
    ncols, nrows = lat.dimensions(5.0, 5.0)
    bx, by, bz = lat.final_box(ncols, nrows, 10.0)
    # final box is close to but not equal to the requested 5.0 x 5.0
    assert abs(bx - 5.0) < 0.6 and bx != 5.0
    assert abs(by - 5.0) < 0.6 and by != 5.0
    assert math.isclose(bx, ncols * lat.colsep)
    assert math.isclose(by, nrows * lat.rowsep)


def test_tilt_preserves_shape():
    # tilting must not change internal distances (rigid rotation)
    coords = np.array([[0, 0, 0.0], [0, 0, 1.0], [0.5, 0, 0.5]])
    out = orient.apply_tilt(coords, alpha=28, beta=53)
    d_in = np.linalg.norm(coords[0] - coords[1])
    d_out = np.linalg.norm(out[0] - out[1])
    assert math.isclose(d_in, d_out, abs_tol=1e-9)


def test_check_oriented_rejects_horizontal():
    # a chain lying along x should fail the sanity check
    coords = np.array([[0, 0, 0.0], [1.0, 0, 0], [2.0, 0, 0]])
    try:
        orient.check_oriented(coords)
        assert False, "expected mis-orientation error"
    except ValueError:
        pass


def test_canonicalize_points_anchor_to_head_up():
    coords = np.array([[0, 0, 0.0], [1.0, 0, 0]])  # anchor->head along +x
    out = orient.canonicalize(coords, anchor_idx=0, head_idx=1)
    v = out[1] - out[0]
    v = v / np.linalg.norm(v)
    assert v[2] > 0.99  # now points along +z


def test_site_density_and_sites_for():
    lat = Lattice.rounded()                      # colsep=0.499, rowsep=0.432
    assert math.isclose(lat.site_density(), 1.0/(0.499*0.432))
    sites = list(lat.sites_for(4, 6))            # explicit grid
    assert len(sites) == 24
    # row 1 is x-offset (hex stagger), row 0 is not
    r0 = [s for s in sites if s[0] == 0][0]
    r1 = [s for s in sites if s[0] == 1][0]
    assert math.isclose(r0[2], 0.0)
    assert math.isclose(r1[2], lat.offset)
