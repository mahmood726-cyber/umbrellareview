# E156 Protocol — `UmbrellaReview`

This repository is the source code and dashboard backing an E156 micro-paper on the [E156 Student Board](https://mahmood726-cyber.github.io/e156/students.html).

---

## `[164]` Umbrella Review Evidence Classifier with Credibility Grading

**Type:** methods  |  ESTIMAND: Evidence conviction class  
**Data:** Summary-level data from multiple meta-analyses

### 156-word body

How can researchers conducting umbrella reviews classify evidence strength across multiple meta-analyses using standardized credibility criteria? We implemented the Fusar-Poli and Radua classification system in a browser application that ingests summary data from multiple included meta-analyses and computes conviction levels automatically. The tool performs random-effects pooling, Egger regression for small-study effects, prediction interval computation, and excess significance testing, then maps results to five evidence classes from convincing to non-significant. Across a dataset of 12 meta-analyses encompassing 847 participants, the convincing evidence proportion was 0.17 (95% CI 0.04 to 0.35), while 4 associations were downgraded to weak. Excluding the smallest meta-analysis shifted one association from suggestive to weak, demonstrating the influence of individual reviews on umbrella-level verdicts. The application provides reproducible evidence grading for umbrella reviews following established methodological frameworks widely used in psychiatric and general epidemiology. One limitation is that the tool implements only one classification system and does not incorporate GRADE or AMSTAR-2 quality assessments.

### Submission metadata

```
Corresponding author: Mahmood Ahmad <mahmood.ahmad2@nhs.net>
ORCID: 0000-0001-9107-3704
Affiliation: Tahir Heart Institute, Rabwah, Pakistan

Links:
  Code:      https://github.com/mahmood726-cyber/UmbrellaReview
  Protocol:  https://github.com/mahmood726-cyber/UmbrellaReview/blob/main/E156-PROTOCOL.md
  Dashboard: https://mahmood726-cyber.github.io/umbrellareview/

References (topic pack: GRADE / certainty rating):
  1. Guyatt GH, Oxman AD, Vist GE, et al. 2008. GRADE: an emerging consensus on rating quality of evidence and strength of recommendations. BMJ. 336(7650):924-926. doi:10.1136/bmj.39489.470347.AD
  2. Schünemann HJ, Cuello C, Akl EA, et al. 2019. GRADE guidelines: 18. How ROBINS-I and other tools to assess risk of bias in nonrandomized studies should be used to rate the certainty of a body of evidence. J Clin Epidemiol. 111:105-114. doi:10.1016/j.jclinepi.2018.01.012

Data availability: No patient-level data used. Analysis derived exclusively
  from publicly available aggregate records. All source identifiers are in
  the protocol document linked above.

Ethics: Not required. Study uses only publicly available aggregate data; no
  human participants; no patient-identifiable information; no individual-
  participant data. No institutional review board approval sought or required
  under standard research-ethics guidelines for secondary methodological
  research on published literature.

Funding: None.

Competing interests: MA serves on the editorial board of Synthēsis (the
  target journal); MA had no role in editorial decisions on this
  manuscript, which was handled by an independent editor of the journal.

Author contributions (CRediT):
  [STUDENT REWRITER, first author] — Writing – original draft, Writing –
    review & editing, Validation.
  [SUPERVISING FACULTY, last/senior author] — Supervision, Validation,
    Writing – review & editing.
  Mahmood Ahmad (middle author, NOT first or last) — Conceptualization,
    Methodology, Software, Data curation, Formal analysis, Resources.

AI disclosure: Computational tooling (including AI-assisted coding via
  Claude Code [Anthropic]) was used to develop analysis scripts and assist
  with data extraction. The final manuscript was human-written, reviewed,
  and approved by the author; the submitted text is not AI-generated. All
  quantitative claims were verified against source data; cross-validation
  was performed where applicable. The author retains full responsibility for
  the final content.

Preprint: Not preprinted.

Reporting checklist: PRISMA 2020 (methods-paper variant — reports on review corpus).

Target journal: ◆ Synthēsis (https://www.synthesis-medicine.org/index.php/journal)
  Section: Methods Note — submit the 156-word E156 body verbatim as the main text.
  The journal caps main text at ≤400 words; E156's 156-word, 7-sentence
  contract sits well inside that ceiling. Do NOT pad to 400 — the
  micro-paper length is the point of the format.

Manuscript license: CC-BY-4.0.
Code license: MIT.

SUBMITTED: [ ]
```


---

_Auto-generated from the workbook by `C:/E156/scripts/create_missing_protocols.py`. If something is wrong, edit `rewrite-workbook.txt` and re-run the script — it will overwrite this file via the GitHub API._