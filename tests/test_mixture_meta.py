"""Tests for finite mixture meta-meta-analysis."""

import math
import pytest
from umbrella.models import ReviewInput
from umbrella.mixture_meta import compute_mixture_meta


def _make_reviews(thetas, ses=None):
    """Helper to build ReviewInput list from theta values."""
    if ses is None:
        ses = [0.15] * len(thetas)
    return [
        ReviewInput(
            review_id=f"R{i+1}", theta=t, ci_lo=t - 1.96 * s, ci_hi=t + 1.96 * s,
            se=s, k=5, study_ids=[f"s{i}_{j}" for j in range(5)], measure="logOR",
        )
        for i, (t, s) in enumerate(zip(thetas, ses))
    ]


def test_components_ge_one(statin_reviews):
    """Number of selected components must be >= 1."""
    result = compute_mixture_meta(statin_reviews)
    assert result["n_components"] >= 1


def test_weights_sum_to_one(statin_reviews):
    """Mixing weights must sum to 1."""
    result = compute_mixture_meta(statin_reviews)
    total = sum(result["weights"])
    assert abs(total - 1.0) < 1e-6, f"Weights sum to {total}"


def test_posterior_probs_per_review(sglt2i_reviews):
    """Each review has posterior probs summing to ~1."""
    result = compute_mixture_meta(sglt2i_reviews)
    for rid, probs in result["posterior_probs"].items():
        total = sum(probs)
        assert abs(total - 1.0) < 1e-6, f"Posterior for {rid} sums to {total}"


def test_bic_finite(statin_reviews):
    """BIC values for all K must be finite."""
    result = compute_mixture_meta(statin_reviews)
    for k, bic in result["bic_by_k"].items():
        assert math.isfinite(bic), f"BIC for K={k} is {bic}"


def test_bimodal_detects_two_components():
    """Two well-separated clusters should favour K>=2."""
    # Cluster 1: theta ~ -0.5, Cluster 2: theta ~ +0.5
    thetas = [-0.6, -0.5, -0.4, -0.55, 0.4, 0.5, 0.6, 0.55]
    reviews = _make_reviews(thetas, ses=[0.05] * 8)
    result = compute_mixture_meta(reviews, max_k=3)
    # Should detect at least 2 components (BIC for K=2 < BIC for K=1)
    assert result["bic_by_k"][2] < result["bic_by_k"][1], \
        f"BIC(K=2)={result['bic_by_k'][2]} not < BIC(K=1)={result['bic_by_k'][1]}"
