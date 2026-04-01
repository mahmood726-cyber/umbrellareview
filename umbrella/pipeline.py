from umbrella.models import UmbrellaVerdict
from umbrella.overlap import compute_overlap
from umbrella.amstar import score_amstar
from umbrella.concordance import compute_concordance
from umbrella.discordance import compute_discordance
from umbrella.verdict import synthesize_recommendation
from umbrella.certifier import compute_input_hash, certify

def run_umbrella(reviews):
    """End-to-end umbrella review analysis."""
    overlap = compute_overlap(reviews)

    amstar_results = []
    for r in reviews:
        if r.amstar_items:
            amstar_results.append(score_amstar(r.review_id, r.amstar_items))
        else:
            from umbrella.models import AMSTARResult
            amstar_results.append(AMSTARResult(review_id=r.review_id, confidence="Low"))

    concordance = compute_concordance(reviews)
    discordance = compute_discordance(reviews, overlap, concordance)
    recommendation = synthesize_recommendation(reviews, overlap, amstar_results, concordance, discordance)

    input_hash = compute_input_hash(reviews)
    cert = certify(reviews, overlap, concordance)

    return UmbrellaVerdict(
        overlap=overlap,
        amstar_results=amstar_results,
        concordance=concordance,
        discordance=discordance,
        recommendation=recommendation,
        input_hash=input_hash,
        certification=cert,
    )
