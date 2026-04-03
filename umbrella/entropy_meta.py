"""Information-Theoretic Meta-Meta-Analysis.

Entropy-based analysis of the systematic review evidence base:
Shannon entropy of review conclusions, mutual information between
quality and conclusion, information gain per review, redundancy,
and Simpson's diversity index.
"""

from __future__ import annotations

import math
from collections import Counter


# ── Quality classification ────────────────────────────────────────

def _amstar_confidence(review) -> str:
    """Derive AMSTAR-2 confidence label from amstar_items dict."""
    items = getattr(review, "amstar_items", None) or {}
    if not items:
        return "Unknown"
    critical = {"1", "2", "4", "7", "9", "11", "13"}
    n_crit_flaw = sum(1 for k, v in items.items() if k in critical and v == "no")
    n_noncrit_weak = sum(1 for k, v in items.items() if k not in critical and v != "yes")
    if n_crit_flaw == 0:
        return "High" if n_noncrit_weak <= 1 else "Moderate"
    if n_crit_flaw == 1:
        return "Low"
    return "CritLow"


def _classify_conclusion(review) -> str:
    """Classify a review as beneficial / harmful / inconclusive from its CI."""
    if review.theta > 0 and review.ci_lo > 0:
        return "beneficial"
    elif review.theta < 0 and review.ci_hi < 0:
        return "harmful"
    return "inconclusive"


# ── Entropy helpers ───────────────────────────────────────────────

def _shannon_entropy(counts: dict) -> float:
    """Shannon entropy H = -sum(p * log2(p)) from a counter dict."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def _joint_entropy(contingency: dict) -> float:
    """H(Q, C) from joint counts."""
    total = sum(contingency.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in contingency.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


# ── Public API ────────────────────────────────────────────────────

def compute_entropy_meta(reviews):
    """Run information-theoretic meta-meta-analysis.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with keys:
        conclusion_entropy           — float (Shannon H in bits)
        normalized_entropy           — float (H / H_max, in [0,1])
        mutual_info_quality_conclusion — float (I(Q;C) in bits)
        conditional_entropy          — float (H(C|Q) in bits)
        info_gain_per_review         — dict review_id -> delta_H
        redundancy                   — float (1 - H/H_max)
        simpson_diversity            — float (1 - sum(p_i^2))
        conclusion_distribution      — dict category -> count
    """
    if not reviews:
        return {
            "conclusion_entropy": 0.0,
            "normalized_entropy": 0.0,
            "mutual_info_quality_conclusion": 0.0,
            "conditional_entropy": 0.0,
            "info_gain_per_review": {},
            "redundancy": 1.0,
            "simpson_diversity": 0.0,
            "conclusion_distribution": {},
        }

    n = len(reviews)

    # ── Classify each review ──
    conclusions = [_classify_conclusion(r) for r in reviews]
    qualities = [_amstar_confidence(r) for r in reviews]

    # ── Conclusion distribution ──
    conc_counts = Counter(conclusions)
    conclusion_distribution = dict(conc_counts)

    # ── Shannon entropy of conclusions ──
    conclusion_entropy = _shannon_entropy(conc_counts)

    # Max entropy: all 3 categories present -> log2(n_categories)
    n_categories = len(set(["beneficial", "harmful", "inconclusive"]))
    h_max = math.log2(n_categories)  # log2(3) ~= 1.585
    normalized_entropy = conclusion_entropy / h_max if h_max > 0 else 0.0

    # ── Mutual information I(Q; C) ──
    # Build contingency table
    quality_counts = Counter(qualities)
    joint_counts = Counter(zip(qualities, conclusions))

    # I(Q;C) = H(Q) + H(C) - H(Q,C)
    h_q = _shannon_entropy(quality_counts)
    h_c = conclusion_entropy
    h_qc = _joint_entropy(joint_counts)
    mutual_info = max(0.0, h_q + h_c - h_qc)

    # ── Conditional entropy H(C|Q) = H(C) - I(Q;C) ──
    conditional_entropy = max(0.0, h_c - mutual_info)

    # ── Information gain per review ──
    info_gain = {}
    for idx, r in enumerate(reviews):
        # Compute entropy without this review
        conc_without = [c for j, c in enumerate(conclusions) if j != idx]
        counts_without = Counter(conc_without)
        h_without = _shannon_entropy(counts_without)
        delta_h = conclusion_entropy - h_without
        info_gain[r.review_id] = round(delta_h, 6)

    # ── Redundancy ──
    redundancy = 1.0 - normalized_entropy

    # ── Simpson's diversity ──
    total = sum(conc_counts.values())
    if total > 0:
        simpson = 1.0 - sum((c / total) ** 2 for c in conc_counts.values())
    else:
        simpson = 0.0

    return {
        "conclusion_entropy": round(conclusion_entropy, 6),
        "normalized_entropy": round(normalized_entropy, 6),
        "mutual_info_quality_conclusion": round(mutual_info, 6),
        "conditional_entropy": round(conditional_entropy, 6),
        "info_gain_per_review": info_gain,
        "redundancy": round(redundancy, 6),
        "simpson_diversity": round(simpson, 6),
        "conclusion_distribution": conclusion_distribution,
    }
