import numpy as np
from umbrella.models import DiscordanceResult, DiscordanceFactor
from umbrella.amstar import score_amstar

def _scope_divergence(reviews):
    """Jaccard distance of scope_tags across reviews."""
    if len(reviews) < 2:
        return 0.0
    tags = [set(r.scope_tags) for r in reviews]
    dists = []
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            union = tags[i] | tags[j]
            inter = tags[i] & tags[j]
            if len(union) == 0:
                dists.append(0.0)
            else:
                dists.append(1.0 - len(inter) / len(union))
    return float(np.mean(dists)) if dists else 0.0

def _inclusion_divergence(reviews):
    """k ratio: max(k)/min(k). High ratio = different inclusion criteria."""
    ks = [r.k for r in reviews if r.k > 0]
    if len(ks) < 2:
        return 0.0
    ratio = max(ks) / min(ks)
    return min(1.0, (ratio - 1.0) / 3.0)  # normalize: ratio 4.0 -> 1.0

def _quality_divergence(reviews):
    """Range of AMSTAR-2 confidence across reviews."""
    levels = {"High": 4, "Moderate": 3, "Low": 2, "Critically Low": 1}
    scores = []
    for r in reviews:
        if r.amstar_items:
            res = score_amstar(r.review_id, r.amstar_items)
            scores.append(levels.get(res.confidence, 2))
        else:
            scores.append(2)
    if len(scores) < 2:
        return 0.0
    return (max(scores) - min(scores)) / 3.0  # normalize: range 3 -> 1.0

def _method_divergence(reviews):
    """Do reviews use different effect measures?"""
    measures = set(r.measure for r in reviews)
    if len(measures) <= 1:
        return 0.0
    return min(1.0, (len(measures) - 1) / 2.0)

def compute_discordance(reviews, overlap, concordance):
    """Classify discordance level and decompose contributing factors."""
    # Classification
    dir_agree = concordance.direction_agreement
    ci_overlap = concordance.ci_overlap_fraction
    i2 = concordance.meta_meta_i2

    # Check for contradictory: opposite significant results
    sig_positive = sum(1 for r in reviews if r.ci_lo > 0)
    sig_negative = sum(1 for r in reviews if r.ci_hi < 0)
    contradictory = sig_positive > 0 and sig_negative > 0

    if contradictory:
        overall = "Contradictory"
    elif dir_agree < 1.0 or i2 > 50:
        overall = "Major"
    elif ci_overlap < 1.0:
        overall = "Minor"
    else:
        overall = "Concordant"

    # Factor decomposition
    scope_d = _scope_divergence(reviews)
    inclusion_d = _inclusion_divergence(reviews)
    overlap_d = 1.0 - overlap.cca  # low CCA = different evidence bases
    quality_d = _quality_divergence(reviews)
    method_d = _method_divergence(reviews)

    raw = [
        ("Scope differences", scope_d, "Reviews differ in population, intervention, or outcome definitions"),
        ("Inclusion criteria", inclusion_d, "Reviews include different numbers of studies (k ratio)"),
        ("Evidence base", overlap_d, "Reviews analyse largely different primary studies (low CCA)"),
        ("Quality variation", quality_d, "Reviews vary in AMSTAR-2 confidence rating"),
        ("Statistical methods", method_d, "Reviews use different effect measures or pooling methods"),
    ]

    total = sum(r[1] for r in raw)
    factors = []
    for name, score, desc in sorted(raw, key=lambda x: -x[1]):
        pct = round(score / total * 100, 1) if total > 0 else 0.0
        factors.append(DiscordanceFactor(factor=name, contribution=pct, description=desc))

    return DiscordanceResult(overall_discordance=overall, factors=factors)
