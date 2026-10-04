"""The locus passport: a build-independent description of one GWAS locus.

Identity comes from the sequence anchors (and the credible set). Everything
else - placements, peak, LD block, genes, neighbours - describes the locus as
seen in one panel and one build, and may differ between studies.

Coordinates are 1-based and inclusive.
"""

from __future__ import annotations

import json
from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .ids import normalize_sequence, passport_id

SCHEMA_VERSION = "0.1.0"


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class AnchorRole(str, Enum):
    lead_snp = "lead_snp"
    block_left = "block_left"
    block_right = "block_right"
    gene = "gene"


class Anchor(_Model):
    """A flanking sequence that locates the locus on any assembly by alignment.

    `variant_offset` marks the position the anchor stands for: the lead
    variant, a block edge, or the centre of a gene anchor.
    """

    role: AnchorRole
    sequence: str
    variant_offset: Optional[int] = Field(
        None, ge=0, description="0-based index of the variant inside `sequence` (lead_snp only)"
    )
    unique: Optional[bool] = Field(None, description="maps once to the source genome; None = not checked")
    label: Optional[str] = Field(None, description="e.g. gene ID for a gene anchor")

    @field_validator("sequence")
    @classmethod
    def _check_sequence(cls, v: str) -> str:
        return normalize_sequence(v)

    @model_validator(mode="after")
    def _check_offset(self) -> Anchor:
        if self.variant_offset is None:
            if self.role is AnchorRole.lead_snp:
                raise ValueError("a lead_snp anchor needs variant_offset")
            self.variant_offset = len(self.sequence) // 2
        if self.variant_offset is not None and self.variant_offset >= len(self.sequence):
            raise ValueError("variant_offset is outside the anchor sequence")
        return self


class LocusKind(str, Enum):
    gwas = "gwas"  # a GWAS peak; the primary anchor is the lead SNP
    gene = "gene"  # a cloned / known gene; the primary anchor is in the gene
    qtl = "qtl"  # a mapped QTL interval with a marker as primary anchor


class PlacementFlag(str, Enum):
    missing_anchor = "missing_anchor"
    multi_mapping = "multi_mapping"
    split = "split"
    inverted = "inverted"
    out_of_order = "out_of_order"


class Placement(_Model):
    """Where the locus sits on one genome build. Derived from anchors, never the identity."""

    build: str
    chrom: str
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    strand: str = Field("+", pattern=r"^[+-]$")
    method: str = Field("source", description="'source' = coordinates from the study; 'anchor' = placed by alignment")
    lead_pos: Optional[int] = Field(None, ge=1, description="position of the lead variant (or gene anchor centre) on this build")
    anchors_placed: Optional[int] = Field(None, ge=0)
    anchors_total: Optional[int] = Field(None, ge=0)
    flags: list[PlacementFlag] = []

    @model_validator(mode="after")
    def _check_interval(self) -> Placement:
        if self.end < self.start:
            raise ValueError("placement end is before start")
        return self


class Variant(_Model):
    id: Optional[str] = None
    chrom: str
    pos: int = Field(ge=1)
    ref: str
    alt: str


class Signal(_Model):
    """The statistical peak, as observed in the source study."""

    lead: Variant
    effect_allele: Optional[str] = None
    effect_allele_freq: Optional[float] = Field(None, ge=0, le=1)
    beta: Optional[float] = None
    se: Optional[float] = Field(None, ge=0)
    pvalue: float = Field(gt=0, le=1)
    peak_start: Optional[int] = Field(None, ge=1, description="first position above the threshold")
    peak_end: Optional[int] = Field(None, ge=1)
    n_significant: Optional[int] = Field(None, ge=0)


class CredibleVariant(_Model):
    id: Optional[str] = None
    pos: int = Field(ge=1)
    pip: float = Field(ge=0, le=1, description="posterior inclusion probability")


class CredibleSet(_Model):
    method: str = Field(description="e.g. 'wakefield_abf', 'susie'")
    coverage: float = Field(0.95, gt=0, le=1)
    variants: list[CredibleVariant]


class LDBlock(_Model):
    """LD around the lead in the source panel. Panel-specific by nature."""

    start: int = Field(ge=1)
    end: int = Field(ge=1)
    method: str = Field(description="e.g. 'r2>=0.6 from lead', 'gabriel'")
    r2_threshold: Optional[float] = Field(None, ge=0, le=1)
    n_haplotypes: Optional[int] = Field(None, ge=1)
    haplotype_freqs: list[float] = []
    inversion_like: Optional[bool] = Field(None, description="long block with near-perfect LD")

    @property
    def length_kb(self) -> float:
        return (self.end - self.start + 1) / 1000


class Imprint(_Model):
    """The locus's surroundings as seen in the GWAS summary statistics of the source study.

    `pos`, `ids`, `ref`, `alt`, `z` are parallel lists for every tested SNP within
    +- window_bp of the lead (source build). `z` is signed for the alt allele when the
    study reports effects (`z_signed`), otherwise |z| from the p-value.
    `profile` is the max |z| in `bins` equal bins across the window, centred on the
    lead: a fixed-length shape that can be compared between studies.
    `ld_edges` links the top significant SNPs (indices into the lists) with r2 >= ld_edge_r2.
    """

    window_bp: int = Field(ge=1)
    threshold: float
    z_signed: bool
    ids: list[str]
    pos: list[int]
    ref: list[str]
    alt: list[str]
    z: list[float]
    bins: int
    profile: list[float]
    half_max_width_bp: int = Field(ge=0, description="span of SNPs with -log10 p >= half the lead's")
    n_significant: int = Field(ge=0)
    independent_signals: int = Field(ge=1, description="loci of this study within the window, this one included")
    ld_edge_r2: Optional[float] = None
    ld_edges: list[tuple[int, int, float]] = []
    territory_start: Optional[int] = Field(None, description="this locus's own part of the window: halfway to the "
                                           "neighbouring leads of the same study (source build)")
    territory_end: Optional[int] = None


class Gene(_Model):
    id: str
    name: Optional[str] = None
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    distance_bp: int = Field(ge=0, description="0 when the gene overlaps the lead or credible set")
    in_credible_set: bool = False
    annotation: Optional[str] = None


class Neighbour(_Model):
    passport_id: str
    distance_bp: int = Field(ge=0)
    r2_with_lead: Optional[float] = Field(None, ge=0, le=1)
    independent: Optional[bool] = None


class Trait(_Model):
    name: str
    ontology_id: Optional[str] = Field(None, description="e.g. Planteome TO:0000207")


class Source(_Model):
    study: str
    panel: Optional[str] = None
    n_samples: Optional[int] = Field(None, ge=1)
    gwas_method: Optional[str] = None
    doi: Optional[str] = None


class LocusPassport(_Model):
    schema_version: str = SCHEMA_VERSION
    passport_id: str = ""
    kind: LocusKind = LocusKind.gwas
    species: str
    trait: Trait
    traits_other: list[Trait] = Field([], description="further traits linked to this locus (pleiotropy)")
    source: Source
    anchors: list[Anchor]
    placements: list[Placement] = []
    signal: Optional[Signal] = None
    credible_set: Optional[CredibleSet] = None
    ld_block: Optional[LDBlock] = None
    imprint: Optional[Imprint] = None
    genes: list[Gene] = []
    neighbours: list[Neighbour] = []
    notes: Optional[str] = None

    @property
    def primary_role(self) -> AnchorRole:
        return AnchorRole.gene if self.kind is LocusKind.gene else AnchorRole.lead_snp

    @model_validator(mode="after")
    def _assign_id(self) -> LocusPassport:
        role = self.primary_role
        primary = [a for a in self.anchors if a.role is role]
        if len(primary) != 1:
            raise ValueError(f"a {self.kind.value} passport needs exactly one {role.value} anchor")
        if self.kind is LocusKind.gwas and self.signal is None:
            raise ValueError("a gwas passport needs a signal")
        expected = passport_id(self.species, primary[0].sequence)
        if not self.passport_id:
            self.passport_id = expected
        elif self.passport_id != expected:
            raise ValueError(f"passport_id {self.passport_id} does not match its lead anchor ({expected})")
        return self

    @property
    def primary_anchor(self) -> Anchor:
        return next(a for a in self.anchors if a.role is self.primary_role)

    @property
    def lead_anchor(self) -> Anchor:
        return self.primary_anchor

    def placement(self, build: str) -> Optional[Placement]:
        return next((p for p in self.placements if p.build == build), None)

    def to_json(self, path: str | Path | None = None, indent: int = 2) -> str:
        text = self.model_dump_json(indent=indent, exclude_none=True)
        if path is not None:
            Path(path).write_text(text + "\n")
        return text

    @classmethod
    def from_json(cls, path: str | Path) -> LocusPassport:
        return cls.model_validate_json(Path(path).read_text())


def json_schema() -> dict:
    return LocusPassport.model_json_schema()


def write_json_schema(path: str | Path) -> None:
    Path(path).write_text(json.dumps(json_schema(), indent=2) + "\n")


def read_passports(path: str | Path) -> list[LocusPassport]:
    """Read a .jsonl collection (one passport per line) or a single .json passport."""
    path = Path(path)
    if path.suffix == ".jsonl":
        return [LocusPassport.model_validate_json(line) for line in path.read_text().splitlines() if line.strip()]
    return [LocusPassport.from_json(path)]


def write_passports(path: str | Path, passports: list[LocusPassport]) -> None:
    with open(path, "w") as fh:
        for p in passports:
            fh.write(p.model_dump_json(exclude_none=True) + "\n")
