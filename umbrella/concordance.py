import numpy as np
from umbrella.models import ConcordanceResult
from umbrella.amstar import score_amstar

AMSTAR_WEIGHTS = {"High": 1.0, "Moderate": 0.75, "Low": 0.5, "Critically Low": 0.25}

def compute_concordance(reviews):
    """Compute effect concordance across reviews."""
    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)  # fallback

    # Direction agreement
    majority_sign = np.sign(np.median(thetas))
    if majority_sign == 0:
        majority_sign = -1.0
    direction_agreement = float(np.mean(np.sign(thetas) == majority_sign))

    # All significant?
    ci_los = np.array([r.ci_lo for r in reviews])
    ci_his = np.array([r.ci_hi for r in reviews])
    all_sig = bool(np.all((ci_los > 0) | (ci_his < 0)))

    # CI overlap fraction (pairwise)
    n = len(reviews)
    overlap_count = 0
    pair_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            lo = max(ci_los[i], ci_los[j])
            hi = min(ci_his[i], ci_his[j])
            if lo < hi:
                overlap_count += 1
            pair_count += 1
    ci_overlap = overlap_count / pair_count if pair_count > 0 else 1.0

    # Meta-meta-analysis (inverse-variance)
    w = 1.0 / (ses ** 2)
    meta_theta = float(np.sum(w * thetas) / np.sum(w))
    meta_se = float(1.0 / np.sqrt(np.sum(w)))
    q = float(np.sum(w * (thetas - meta_theta) ** 2))
    i2 = max(0.0, (q - (n - 1)) / q * 100) if q > 0 and n > 1 else 0.0

    # Quality-weighted
    amstar_w = []
    for r in reviews:
        if r.amstar_items:
            res = score_amstar(r.review_id, r.amstar_items)
            amstar_w.append(AMSTAR_WEIGHTS.get(res.confidence, 0.5))
        else:
            amstar_w.append(0.5)
    amstar_w = np.array(amstar_w)
    w_adj = amstar_w / (ses ** 2)
    qw_theta = float(np.sum(w_adj * thetas) / np.sum(w_adj))
    qw_se = float(1.0 / np.sqrt(np.sum(w_adj)))
    qw_ci = (qw_theta - 1.96 * qw_se, qw_theta + 1.96 * qw_se)

    return ConcordanceResult(
        direction_agreement=round(direction_agreement, 3),
        all_significant=all_sig,
        ci_overlap_fraction=round(ci_overlap, 3),
        meta_meta_theta=round(meta_theta, 4),
        meta_meta_se=round(meta_se, 4),
        meta_meta_i2=round(i2, 1),
        range_theta=(round(float(thetas.min()), 4), round(float(thetas.max()), 4)),
        quality_weighted_theta=round(qw_theta, 4),
        quality_weighted_ci=(round(qw_ci[0], 4), round(qw_ci[1], 4)),
    )
