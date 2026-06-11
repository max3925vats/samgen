from tests._synthetic import linear_strand
from samgen.geometry import generate_geometry

def _tile(mol, M, tmp_path):
    cfg = {"lattice": {"rounded": True, "tilt_alpha": 0, "tilt_beta": 0},
           "box": {"x": 3.0, "y": 3.0, "z": 8.0},
           "design": {"type": "uniform", "component": "c"}}
    out = str(tmp_path / "s.gro")
    res = generate_geometry(cfg, {"c": mol}, out_gro=out, is_tty=False)
    strands = res.manifest["counts"]["c"]
    assert res.manifest["natoms"] == strands * M

def test_placement_is_size_agnostic(tmp_path):
    for M in (7, 33):
        mol = linear_strand("X", n_chain=M - 1, with_cap=False)  # M atoms total
        _tile(mol, M, tmp_path)
