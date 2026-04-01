
def synthesize_recommendation(reviews, overlap, amstar_results, concordance, discordance):
    """Generate a human-readable recommendation."""
    n = len(reviews)
    dir_pct = round(concordance.direction_agreement * 100)
    theta = concordance.quality_weighted_theta
    ci = concordance.quality_weighted_ci

    # Median AMSTAR confidence
    conf_order = ["Critically Low", "Low", "Moderate", "High"]
    confs = [a.confidence for a in amstar_results]
    conf_ranks = sorted([conf_order.index(c) for c in confs])
    median_conf = conf_order[conf_ranks[len(conf_ranks) // 2]]

    disc = discordance.overall_discordance
    top_factor = discordance.factors[0].factor if discordance.factors else "none identified"

    parts = []
    parts.append(f"{dir_pct}% of {n} reviews agree on effect direction.")
    parts.append(f"Quality-weighted pooled effect: {theta:.3f} ({ci[0]:.3f} to {ci[1]:.3f}).")
    parts.append(f"Study overlap: CCA={overlap.cca:.3f} ({overlap.groove}).")
    parts.append(f"Median AMSTAR-2 confidence: {median_conf}.")
    parts.append(f"Discordance: {disc}.")
    if disc != "Concordant":
        parts.append(f"Primary driver: {top_factor} ({discordance.factors[0].contribution:.0f}%).")

    return " ".join(parts)
