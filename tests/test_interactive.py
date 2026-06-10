"""Interactive anchor resolution: prompt-first / consent-to-guess.

Uses the bundled ch3 example strand (anchor S12, the 12th atom).
"""

import os
import pytest

from samgen.core.molecule import Molecule
from samgen.interactive import resolve_anchor_interactive

ROOT = os.path.join(os.path.dirname(__file__), "..")
INPUTS = os.path.join(ROOT, "configs", "inputs")
CH3_GRO = os.path.join(INPUTS, "ch3.gro")
CH3_ITP = os.path.join(INPUTS, "ch3.itp")

S12_IDX = 11  # S12 is the 12th atom (0-based 11)


def _ch3():
    return Molecule.from_files("CH3", CH3_GRO, CH3_ITP)


def _replies(*answers):
    it = iter(answers)
    return lambda prompt: next(it)


def test_specified_anchor_does_not_prompt():
    def boom(prompt):  # must never be called when anchor is specified
        raise AssertionError("should not prompt")
    res = resolve_anchor_interactive(_ch3(), "S12", False, input_fn=boom, is_tty=True)
    assert res.anchor_idx == S12_IDX


def test_batch_without_consent_raises():
    with pytest.raises(ValueError, match="no anchor specified"):
        resolve_anchor_interactive(_ch3(), None, allow_autodetect=False,
                                   input_fn=_replies(), is_tty=False)


def test_batch_with_consent_autodetects():
    res = resolve_anchor_interactive(_ch3(), None, allow_autodetect=True,
                                     input_fn=_replies(), is_tty=False)
    assert res.cap_carbon_idx is not None  # found S12 + its methyl cap


def test_interactive_typed_name():
    res = resolve_anchor_interactive(_ch3(), None, False,
                                     input_fn=_replies("S12"), is_tty=True)
    assert res.anchor_idx == S12_IDX


def test_interactive_blank_then_consent_then_confirm():
    res = resolve_anchor_interactive(_ch3(), None, False,
                                     input_fn=_replies("", "y", "y"), is_tty=True)
    assert res.anchor_idx == S12_IDX


def test_interactive_declines_autodetect():
    with pytest.raises(ValueError, match="declined"):
        resolve_anchor_interactive(_ch3(), None, False,
                                   input_fn=_replies("", "n"), is_tty=True)


def test_interactive_rejects_guess():
    with pytest.raises(ValueError, match="rejected"):
        resolve_anchor_interactive(_ch3(), None, False,
                                   input_fn=_replies("", "y", "n"), is_tty=True)
