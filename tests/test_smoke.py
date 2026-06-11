import os, subprocess, sys, pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")

def test_uniform_example_runs(tmp_path):
    out = str(tmp_path / "sam.gro")
    r = subprocess.run([sys.executable, "-m", "samgen.cli", "geometry",
                        os.path.join(ROOT, "configs/ch3_onesided.yaml"), "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert os.path.exists(out)
