# UmbrellaEngine

Browser-based tool for conducting umbrella reviews (systematic reviews of systematic reviews). Automates the three hardest parts: study overlap quantification (CCA/GROOVE), quality assessment (AMSTAR-2), and discordance analysis.

## Quick Start

### Browser App (no installation)

Open `app/umbrella.html` in any modern browser. Load a built-in example or enter your own reviews.

### Python Engine

```bash
pip install -e .
```

```python
from umbrella.pipeline import run_umbrella
from umbrella.models import ReviewInput

reviews = [
    ReviewInput(review_id="A", theta=-0.5, ci_lo=-0.8, ci_hi=-0.2,
                k=10, study_ids=["s1","s2","s3"], measure="logOR"),
    ReviewInput(review_id="B", theta=-0.3, ci_lo=-0.6, ci_hi=-0.1,
                k=8, study_ids=["s1","s2","s4"], measure="logOR"),
]
result = run_umbrella(reviews)
print(result.recommendation)
```

## Features

| Feature | Description |
|---------|-------------|
| CCA / GROOVE | Corrected Covered Area with GROOVE classification (Pieper et al. 2014) |
| AMSTAR-2 | 16-item quality assessment with 7 critical domains (Shea et al. 2017) |
| Meta-meta-analysis | Inverse-variance pooling across review-level effects |
| Quality-weighted pooling | AMSTAR-2 confidence as additional weight |
| Discordance decomposition | 5-factor analysis: scope, inclusion, overlap, quality, methods |
| TruthCert | Hash-linked certification bundle for reproducibility |

## Built-in Examples

1. **Statins for Primary CV Prevention** (5 SRs) - Minor discordance
2. **SGLT2i in Heart Failure** (3 SRs) - Concordant
3. **Ivermectin for COVID-19** (4 SRs) - Major discordance, quality-driven

## Browser App Tabs

1. **Input** - Add reviews, load examples, CSV/JSON import
2. **Overlap** - CCA heatmap, study frequency bar chart
3. **Quality** - AMSTAR-2 traffic light grid
4. **Effects** - Forest-of-forests with quality-weighted summary
5. **Discordance** - Waterfall chart of contributing factors
6. **Report & Certify** - Structured report, methods paragraph, TruthCert JSON

## Tests

```bash
python -m pytest tests/ -v
```

## Validation

- CCA formula validated against Pieper et al. 2014
- AMSTAR-2 scoring validated against Shea et al. 2017 decision rules
- Meta-meta-analysis uses standard inverse-variance fixed-effect model

## Citation

Mahmood Ahmad. UmbrellaEngine: A browser-based tool for umbrella reviews. 2026.

## License

MIT
