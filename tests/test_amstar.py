import pytest
from umbrella.amstar import score_amstar, CRITICAL_ITEMS

def _all_yes():
    return {str(i): "yes" for i in range(1, 17)}

def _with_flaws(base, flaws):
    d = dict(base)
    d.update(flaws)
    return d

def test_all_yes_is_high():
    result = score_amstar("R1", _all_yes())
    assert result.confidence == "High"
    assert result.n_critical_flaw == 0

def test_one_critical_no_is_low():
    items = _with_flaws(_all_yes(), {"2": "no"})  # item 2 is critical
    result = score_amstar("R1", items)
    assert result.confidence == "Low"
    assert result.n_critical_flaw == 1

def test_two_critical_no_is_critically_low():
    items = _with_flaws(_all_yes(), {"2": "no", "4": "no"})
    result = score_amstar("R1", items)
    assert result.confidence == "Critically Low"
    assert result.n_critical_flaw == 2

def test_no_critical_two_noncritical_is_moderate():
    items = _with_flaws(_all_yes(), {"3": "no", "5": "no"})
    result = score_amstar("R1", items)
    assert result.confidence == "Moderate"

def test_critical_items_count():
    assert len(CRITICAL_ITEMS) == 7
