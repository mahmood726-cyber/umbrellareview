from umbrella.models import OverlapResult

def _classify_groove(cca):
    if cca < 0.05:
        return "Slight"
    if cca < 0.10:
        return "Moderate"
    if cca < 0.15:
        return "High"
    return "Very High"

def compute_overlap(reviews):
    """Compute CCA, GROOVE, overlap matrix, and study frequency."""
    c = len(reviews)
    if c < 2:
        return OverlapResult(cca=0.0, groove="Slight", n_total_citations=0,
                             n_unique_studies=0, n_reviews=c)

    # Study sets per review
    sets = [set(r.study_ids) for r in reviews]
    all_studies = set()
    for s in sets:
        all_studies |= s

    r = len(all_studies)  # unique studies
    n_total = sum(len(s) for s in sets)  # total citations

    # CCA
    denominator = r * c - r
    if denominator <= 0:
        cca = 0.0
    else:
        cca = (n_total - r) / denominator

    cca = max(0.0, min(1.0, cca))

    # Pairwise overlap matrix
    matrix = [[0] * c for _ in range(c)]
    for i in range(c):
        for j in range(c):
            matrix[i][j] = len(sets[i] & sets[j])

    # Study frequency
    freq = {}
    for sid in all_studies:
        count = sum(1 for s in sets if sid in s)
        freq[sid] = count

    return OverlapResult(
        cca=round(cca, 4),
        groove=_classify_groove(cca),
        n_total_citations=n_total,
        n_unique_studies=r,
        n_reviews=c,
        overlap_matrix=matrix,
        study_frequency=freq,
    )
