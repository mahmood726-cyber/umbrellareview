from dataclasses import dataclass, field
import math

@dataclass
class ReviewInput:
    review_id: str
    theta: float
    ci_lo: float
    ci_hi: float
    se: float = 0.0
    k: int = 0
    study_ids: list[str] = field(default_factory=list)
    measure: str = "OR"
    amstar_items: dict[str, str] = field(default_factory=dict)
    scope_tags: list[str] = field(default_factory=list)
    year: int = 2024
    label: str = ""

    def __post_init__(self):
        if self.se <= 0 and self.ci_lo != 0 and self.ci_hi != 0:
            self.se = (self.ci_hi - self.ci_lo) / (2 * 1.96)
        if not self.label:
            self.label = self.review_id

@dataclass
class OverlapResult:
    cca: float
    groove: str
    n_total_citations: int
    n_unique_studies: int
    n_reviews: int
    overlap_matrix: list[list[int]] = field(default_factory=list)
    study_frequency: dict[str, int] = field(default_factory=dict)

@dataclass
class AMSTARResult:
    review_id: str
    item_scores: dict[str, str] = field(default_factory=dict)
    n_yes: int = 0
    n_critical_yes: int = 0
    n_critical_flaw: int = 0
    confidence: str = ""

@dataclass
class ConcordanceResult:
    direction_agreement: float = 0.0
    all_significant: bool = False
    ci_overlap_fraction: float = 0.0
    meta_meta_theta: float = 0.0
    meta_meta_se: float = 0.0
    meta_meta_i2: float = 0.0
    range_theta: tuple[float, float] = (0.0, 0.0)
    quality_weighted_theta: float = 0.0
    quality_weighted_ci: tuple[float, float] = (0.0, 0.0)

@dataclass
class DiscordanceFactor:
    factor: str
    contribution: float
    description: str

@dataclass
class DiscordanceResult:
    overall_discordance: str = ""
    factors: list[DiscordanceFactor] = field(default_factory=list)

@dataclass
class UmbrellaVerdict:
    overlap: OverlapResult = field(default_factory=lambda: OverlapResult(0, "", 0, 0, 0))
    amstar_results: list[AMSTARResult] = field(default_factory=list)
    concordance: ConcordanceResult = field(default_factory=ConcordanceResult)
    discordance: DiscordanceResult = field(default_factory=DiscordanceResult)
    recommendation: str = ""
    input_hash: str = ""
    certification: str = ""
