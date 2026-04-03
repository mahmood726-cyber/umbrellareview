import pytest
from umbrella.network_meta_meta import network_meta_meta


def test_n_interventions_at_least_two(statin_reviews):
    """Network must identify at least 2 interventions (treatment vs control)."""
    result = network_meta_meta(statin_reviews)
    assert result["n_interventions"] >= 2


def test_p_scores_in_valid_range(sglt2i_reviews):
    """P-scores must be in [0, 1] for each intervention."""
    result = network_meta_meta(sglt2i_reviews)
    for intv, score in result["p_scores"].items():
        assert 0.0 <= score <= 1.0, f"{intv} P-score {score} out of range"


def test_comparisons_nonempty(ivermectin_reviews):
    """Direct comparisons list must not be empty."""
    result = network_meta_meta(ivermectin_reviews)
    assert len(result["comparisons"]) > 0


def test_network_connected_is_bool(statin_reviews):
    """network_connected must be a boolean."""
    result = network_meta_meta(statin_reviews)
    assert isinstance(result["network_connected"], bool)


def test_comparisons_have_required_keys(sglt2i_reviews):
    """Each comparison must have a, b, theta, se, source."""
    result = network_meta_meta(sglt2i_reviews)
    required = {"a", "b", "theta", "se", "source"}
    for comp in result["comparisons"]:
        assert required.issubset(comp.keys()), f"Missing keys in {comp}"
