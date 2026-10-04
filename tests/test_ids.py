import pytest

from plocust.ids import canonical_sequence, passport_id, reverse_complement, species_code


def test_reverse_complement():
    assert reverse_complement("AACGTN") == "NACGTT"


def test_canonical_is_strand_independent():
    seq = "GATTACAGGC"
    assert canonical_sequence(seq) == canonical_sequence(reverse_complement(seq))


def test_passport_id_is_strand_and_case_independent():
    seq = "ACGTTGCAAGGCTTAC"
    a = passport_id("Oryza sativa", seq)
    assert a == passport_id("Oryza sativa", reverse_complement(seq))
    assert a == passport_id("Oryza sativa", seq.lower())
    assert a.startswith("OsLP.")


def test_passport_id_depends_on_species_and_sequence():
    seq = "ACGTTGCAAGGCTTAC"
    assert passport_id("Oryza sativa", seq) != passport_id("Oryza glaberrima", seq)
    assert passport_id("Oryza sativa", seq) != passport_id("Oryza sativa", seq + "A")


def test_species_code():
    assert species_code("Oryza sativa") == "Os"
    with pytest.raises(ValueError):
        species_code("rice")


def test_rejects_invalid_bases():
    with pytest.raises(ValueError):
        canonical_sequence("ACGU")
