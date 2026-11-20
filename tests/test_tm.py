"""Tests for nearest-neighbor melting temperature."""

from pathlib import Path

import pytest
from Bio.SeqUtils.MeltingTemp import Tm_NN

from dna_sequence_analyzer.gc import gc_percent
from dna_sequence_analyzer.tm import (
    OLIGO_GUIDANCE_LENGTH,
    MeltingConditions,
    analyze_melting,
    melting_temperature,
)

BIOPYTHON_EXAMPLE = "CGTTCCAAAGATGTGGGCATGAGCTTAC"


def test_default_matches_biopython_tm_nn() -> None:
    # DNA_NN3, 50 mM Na+, 25 nM of each strand, salt correction 5.
    assert melting_temperature(BIOPYTHON_EXAMPLE) == pytest.approx(Tm_NN(BIOPYTHON_EXAMPLE))


def test_cations_and_oligo_concentration_raise_tm() -> None:
    # Phosphates on the two backbones repel each other. Na+ and Mg2+
    # screen that repulsion, so Tm goes up. dNTPs bind Mg2+ and give
    # some of that screening back. More oligo also raises Tm, because
    # the duplex is favored when the strands are more concentrated.
    oligo = BIOPYTHON_EXAMPLE
    sodium = melting_temperature(oligo, MeltingConditions(na_mM=50))
    more_sodium = melting_temperature(oligo, MeltingConditions(na_mM=200))
    magnesium = melting_temperature(oligo, MeltingConditions(na_mM=50, mg_mM=2))
    magnesium_bound = melting_temperature(
        oligo, MeltingConditions(na_mM=50, mg_mM=2, dntps_mM=1.5)
    )
    dilute = melting_temperature(oligo, MeltingConditions(dnac1_nM=25, dnac2_nM=25))
    concentrated = melting_temperature(oligo, MeltingConditions(dnac1_nM=250, dnac2_nM=250))
    assert more_sodium > sodium
    assert magnesium > sodium
    assert sodium < magnesium_bound < magnesium
    assert concentrated > dilute


def test_self_complementary_matches_biopython() -> None:
    conditions = MeltingConditions(selfcomp=True)
    assert melting_temperature("ATGCAT", conditions) == pytest.approx(
        Tm_NN("ATGCAT", selfcomp=True)
    )


def test_custom_conditions_are_forwarded() -> None:
    conditions = MeltingConditions(na_mM=100, mg_mM=1.5, dnac1_nM=200, dnac2_nM=0, saltcorr=5)
    expected = Tm_NN(BIOPYTHON_EXAMPLE, Na=100, Mg=1.5, dnac1=200, dnac2=0, saltcorr=5)
    assert melting_temperature(BIOPYTHON_EXAMPLE, conditions) == pytest.approx(expected)


def test_same_gc_content_can_have_different_tm() -> None:
    # Both oligos are 8 nt and 50% GC. The stacks are not the same,
    # so a nearest-neighbor Tm must not collapse them to one number.
    alternating = "ATATATAT"
    blocked = "AAAATTTT"
    assert gc_percent(alternating) == pytest.approx(gc_percent(blocked))
    assert melting_temperature(alternating) != pytest.approx(melting_temperature(blocked))


def test_case_whitespace_and_uracil() -> None:
    plain = melting_temperature("ATGCATGC")
    assert melting_temperature("atg catgc") == pytest.approx(plain)
    assert melting_temperature("AUGCAUGC") == pytest.approx(plain)


def test_ambiguous_bases_are_not_stripped() -> None:
    with pytest.raises(ValueError, match="not removed"):
        melting_temperature("ATGCN")
    with pytest.raises(ValueError, match="not removed"):
        melting_temperature("AT-GC")


def test_short_or_empty_sequence_is_rejected() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        melting_temperature("A")
    with pytest.raises(ValueError, match="at least 2"):
        melting_temperature("   ")


def test_conditions_reject_impossible_concentrations() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        MeltingConditions(na_mM=-1)
    with pytest.raises(ValueError, match="saltcorr"):
        MeltingConditions(saltcorr=8)
    with pytest.raises(ValueError, match="dnac1_nM"):
        MeltingConditions(dnac1_nM=10, dnac2_nM=25)


def test_analyze_melting_keeps_a_bad_record(tmp_path: Path) -> None:
    fasta = tmp_path / "oligos.fasta"
    fasta.write_text(">good\nATGCATGC\n>unknown\nATGCN\n>tiny\nA\n", encoding="utf-8")

    records = analyze_melting(fasta)

    assert [item.id for item in records] == ["good", "unknown", "tiny"]
    assert records[0].tm_celsius == pytest.approx(melting_temperature("ATGCATGC"))
    assert records[0].length == 8
    assert records[0].sequence_is_long is False
    assert records[1].tm_celsius is None
    assert records[1].length == 5
    assert records[2].tm_celsius is None
    assert records[2].length == 1


def test_long_sequence_is_flagged_but_still_calculated(tmp_path: Path) -> None:
    sequence = "ATGC" * (OLIGO_GUIDANCE_LENGTH // 4 + 1)
    assert len(sequence) > OLIGO_GUIDANCE_LENGTH
    fasta = tmp_path / "long.fasta"
    fasta.write_text(f">mrna\n{sequence}\n", encoding="utf-8")

    records = analyze_melting(fasta)

    assert records[0].sequence_is_long is True
    assert records[0].tm_celsius == pytest.approx(Tm_NN(sequence))
