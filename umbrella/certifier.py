import hashlib
import json

def compute_input_hash(reviews):
    data = [{"id": r.review_id, "theta": r.theta, "k": r.k} for r in reviews]
    raw = json.dumps(data, sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]

def certify(reviews, overlap, concordance):
    if len(reviews) < 2:
        return "REJECT"
    if overlap.n_unique_studies < 3:
        return "WARN"
    return "PASS"
