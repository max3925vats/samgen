"""Interactive anchor resolution on a synthetic strand (anchor S = atom 1)."""

import pytest

from samgen.interactive import resolve_anchor_interactive
from tests._synthetic import linear_strand

S_IDX = 0  # the synthetic builder puts the sulfur first


def _strand():
    return linear_strand("LIG", n_chain=11, with_cap=True)


def _replies(*answers):
    it = iter(answers)
    return lambda prompt: next(it)


def test_specified_anchor_does_not_prompt():
    def boom(prompt):
        raise AssertionError("should not prompt")
    res = resolve_anchor_interactive(_strand(), "S1", False, input_fn=boom, is_tty=True)
    assert res.anchor_idx == S_IDX


def test_batch_without_consent_raises():
    with pytest.raises(ValueError, match="no anchor specified"):
        resolve_anchor_interactive(_strand(), None, allow_autodetect=False,
                                   input_fn=_replies(), is_tty=False)


def test_batch_with_consent_autodetects():
    res = resolve_anchor_interactive(_strand(), None, allow_autodetect=True,
                                     input_fn=_replies(), is_tty=False)
    assert res.cap_carbon_idx is not None


def test_interactive_blank_then_consent_then_confirm():
    res = resolve_anchor_interactive(_strand(), None, False,
                                     input_fn=_replies("", "y", "y"), is_tty=True)
    assert res.anchor_idx == S_IDX


def test_interactive_declines_autodetect():
    with pytest.raises(ValueError, match="declined"):
        resolve_anchor_interactive(_strand(), None, False,
                                   input_fn=_replies("", "n"), is_tty=True)


def test_interactive_rejects_guess():
    with pytest.raises(ValueError, match="rejected"):
        resolve_anchor_interactive(_strand(), None, False,
                                   input_fn=_replies("", "y", "n"), is_tty=True)
