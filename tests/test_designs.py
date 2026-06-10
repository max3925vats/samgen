from samgen.design.fraction import Fraction


def test_fraction_is_deterministic_by_seed():
    f1 = Fraction(base="b", ligand="L", fraction=0.3, seed=7)
    f2 = Fraction(base="b", ligand="L", fraction=0.3, seed=7)
    labels1 = [f1.label(r, c) for r in range(20) for c in range(20)]
    labels2 = [f2.label(r, c) for r in range(20) for c in range(20)]
    assert labels1 == labels2                       # reproducible
    frac = labels1.count("L") / len(labels1)
    assert 0.2 < frac < 0.4                          # ~target fraction
