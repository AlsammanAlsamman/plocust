import json

import pytest
from pydantic import ValidationError

from plocust import LocusPassport
from plocust.cli import main
from plocust.ids import reverse_complement

# Toy sequences only; not real rice sequence.
LEAD_SEQ = "ACGTTGCAAGGCTTACGATCGATTACAGGCATCG"


def make(**overrides) -> dict:
    d = {
        "species": "Oryza sativa",
        "trait": {"name": "plant height", "ontology_id": "TO:0000207"},
        "source": {"study": "toy", "panel": "toy panel", "n_samples": 400},
        "anchors": [
            {"role": "lead_snp", "sequence": LEAD_SEQ, "variant_offset": 17},
            {"role": "block_left", "sequence": "GGGCCCAAATTT"},
        ],
        "placements": [{"build": "IRGSP-1.0", "chrom": "chr01", "start": 1000, "end": 25000}],
        "signal": {
            "lead": {"id": "snp1", "chrom": "chr01", "pos": 12000, "ref": "A", "alt": "G"},
            "pvalue": 1e-12,
            "beta": -0.8,
        },
        "ld_block": {"start": 1000, "end": 25000, "method": "r2>=0.6 from lead", "r2_threshold": 0.6},
    }
    d.update(overrides)
    return d


def test_id_assigned_from_lead_anchor():
    p = LocusPassport(**make())
    assert p.passport_id.startswith("OsLP.")
    assert p.lead_anchor.variant_offset == 17
    assert p.ld_block.length_kb == 24.001


def test_id_survives_new_build_and_strand():
    p1 = LocusPassport(**make())
    d = make(placements=[{"build": "MH63", "chrom": "chr01", "start": 900, "end": 24000, "method": "anchor"}])
    d["anchors"][0] = {"role": "lead_snp", "sequence": reverse_complement(LEAD_SEQ), "variant_offset": 16}
    assert LocusPassport(**d).passport_id == p1.passport_id


def test_json_round_trip(tmp_path):
    p = LocusPassport(**make())
    path = tmp_path / "p.json"
    p.to_json(path)
    assert LocusPassport.from_json(path) == p


def test_wrong_id_rejected():
    with pytest.raises(ValidationError, match="does not match"):
        LocusPassport(**make(passport_id="OsLP.wrong"))


@pytest.mark.parametrize(
    "anchors",
    [
        [{"role": "block_left", "sequence": "ACGT"}],  # no lead
        [{"role": "lead_snp", "sequence": "ACGT"}],  # no offset
        [{"role": "lead_snp", "sequence": "ACGT", "variant_offset": 4}],  # offset out of range
    ],
)
def test_bad_anchors_rejected(anchors):
    with pytest.raises(ValidationError):
        LocusPassport(**make(anchors=anchors))


def test_unknown_field_rejected():
    with pytest.raises(ValidationError):
        LocusPassport(**make(colour="red"))


def test_cli_validate(tmp_path, capsys):
    good = tmp_path / "good.json"
    LocusPassport(**make()).to_json(good)
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(make(anchors=[])))
    assert main(["validate", str(good)]) == 0
    assert main(["validate", str(good), str(bad)]) == 1
    assert "FAIL" in capsys.readouterr().err
