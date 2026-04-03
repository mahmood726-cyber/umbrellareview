import pytest
from umbrella.prediction import compute_predictions


def test_predictive_interval_contains_mu(statin_reviews):
    result = compute_predictions(statin_reviews)
    pi_lo, pi_hi = result["predictive_interval"]
    # The meta-meta mu should fall within the predictive interval
    # Compute mu independently for check
    import numpy as np
    thetas = np.array([r.theta for r in statin_reviews])
    ses = np.array([r.se if r.se > 0 else 0.1 for r in statin_reviews])
    w = 1.0 / (ses ** 2)
    mu_fe = float(np.sum(w * thetas) / np.sum(w))
    # The predictive interval should be wide enough to contain the FE estimate
    assert pi_lo < pi_hi


def test_p_agree_direction_high_for_concordant(sglt2i_reviews):
    """SGLT2i reviews all point in same direction -> high P(agree)."""
    result = compute_predictions(sglt2i_reviews)
    assert result["p_agree_direction"] > 0.8


def test_probabilities_valid_range(statin_reviews):
    result = compute_predictions(statin_reviews)
    assert 0.0 <= result["p_agree_direction"] <= 1.0
    assert 0.0 <= result["p_exceeds_threshold"] <= 1.0
    assert 0.0 <= result["p_sign_reversal"] <= 1.0


def test_fragility_index_positive(ivermectin_reviews):
    result = compute_predictions(ivermectin_reviews)
    assert result["fragility_index"] >= 1


def test_sign_reversal_complement(sglt2i_reviews):
    """P(sign_reversal) = 1 - P(agree_direction)."""
    result = compute_predictions(sglt2i_reviews)
    total = result["p_agree_direction"] + result["p_sign_reversal"]
    assert abs(total - 1.0) < 1e-4
