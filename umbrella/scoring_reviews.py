"""Proper Scoring Rules for Review Ranking in Umbrella Reviews."""

import math
import numpy as np
from scipy import stats


def _amstar_weight(review):
    """Compute AMSTAR quality weight in [0, 1]."""
    items = review.amstar_items
    if not items:
        return 0.5  # neutral default when no AMSTAR data
    n_total = len(items)
    if n_total == 0:
        return 0.5
    n_yes = sum(1 for v in items.values() if v.lower() in ("yes", "y"))
    return n_yes / n_total


def compute_scoring(reviews):
    """Compute proper scoring rules and rank reviews.

    Uses the inverse-variance-weighted consensus (meta-meta-analytic mean)
    as the reference value.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with keys: brier_scores, log_scores, crps_scores, composite_scores,
                    ranking, skill_scores, best_review
    """
    n = len(reviews)
    if n < 2:
        raise ValueError("Need at least 2 reviews for scoring")

    thetas = np.array([r.theta for r in reviews], dtype=float)
    ses = np.array([r.se for r in reviews], dtype=float)
    ses = np.where(ses > 0, ses, 0.1)

    # Consensus: inverse-variance weighted mean
    w = 1.0 / (ses ** 2)
    consensus = float(np.sum(w * thetas) / np.sum(w))

    # Range of thetas for Brier normalization
    theta_range = float(np.max(thetas) - np.min(thetas))
    if theta_range < 1e-15:
        theta_range = 1.0  # avoid division by zero when all thetas identical

    brier_scores = {}
    log_scores = {}
    crps_scores = {}
    composite_scores = {}
    skill_scores = {}

    # Reference CRPS: N(0, 100) evaluated at consensus
    sigma_ref = 100.0
    z_ref = (consensus - 0.0) / sigma_ref
    crps_ref = sigma_ref * (
        z_ref * (2 * stats.norm.cdf(z_ref) - 1)
        + 2 * stats.norm.pdf(z_ref)
        - 1.0 / math.sqrt(math.pi)
    )

    for i, r in enumerate(reviews):
        rid = r.review_id
        theta_i = thetas[i]
        se_i = ses[i]

        # --- Brier score ---
        brier = (theta_i - consensus) ** 2 / (theta_range ** 2)
        brier = float(np.clip(brier, 0.0, 1.0))
        brier_scores[rid] = round(brier, 6)

        # --- Log score ---
        # log N(consensus | theta_i, se_i^2)
        log_s = float(stats.norm.logpdf(consensus, loc=theta_i, scale=se_i))
        log_scores[rid] = round(log_s, 6)

        # --- CRPS (closed-form for Normal) ---
        # CRPS = sigma * (z * (2*Phi(z)-1) + 2*phi(z) - 1/sqrt(pi))
        z = (consensus - theta_i) / se_i
        crps = se_i * (
            z * (2 * stats.norm.cdf(z) - 1)
            + 2 * stats.norm.pdf(z)
            - 1.0 / math.sqrt(math.pi)
        )
        crps = float(crps)
        crps_scores[rid] = round(crps, 6)

        # --- Skill score ---
        if crps_ref > 1e-15:
            skill = 1.0 - crps / crps_ref
        else:
            skill = 0.0
        skill_scores[rid] = round(float(skill), 6)

        # --- Composite score ---
        amstar_w = _amstar_weight(r)
        # Lower composite = better (Brier and CRPS are loss functions)
        # Invert amstar_w so that higher quality -> lower composite
        composite = 0.3 * brier + 0.3 * crps + 0.4 * (1 - amstar_w)
        composite_scores[rid] = round(float(composite), 6)

    # Ranking: sort by composite (ascending = best first)
    ranking = sorted(composite_scores.keys(), key=lambda k: composite_scores[k])

    best_review = ranking[0]

    return {
        "brier_scores": brier_scores,
        "log_scores": log_scores,
        "crps_scores": crps_scores,
        "composite_scores": composite_scores,
        "ranking": ranking,
        "skill_scores": skill_scores,
        "best_review": best_review,
    }
