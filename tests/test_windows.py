"""Tests for sliding-window GC, CpG O/E and GC skew."""

from pathlib import Path

import pytest

from dna_sequence_analyzer.gc import gc_extreme_flag, nucleic_acid_kind
from dna_sequence_analyzer.windows import (
    cpg_observed_expected,
    gc_skew,
    gc_windows,
    skew_windows,
)


def test_nucleic_acid_kind_detects_uracil() -> None:
    assert nucleic_acid_kind("ATGC") == "dna"
    assert nucleic_acid_kind("augc") == "rna"
    assert nucleic_acid_kind("ATUC") == "mixed"
    assert nucleic_acid_kind("GGGG") == "unspecified"
    assert nucleic_acid_kind("NNNN") == "unspecified"


def test_extreme_gc_thresholds() -> None:
    assert gc_extreme_flag(24.9) == "low-gc"
    assert gc_extreme_flag(25.0) is None
    assert gc_extreme_flag(65.0) is None
    assert gc_extreme_flag(65.1) == "high-gc"
    assert gc_extreme_flag(None) is None


def test_cpg_ratio_and_skew_are_literal_counts() -> None:
    # CCGG has one CG dinucleotide, two C and two G.
    # O/E = (1 * 4) / (2 * 2) = 1.
    assert cpg_observed_expected("CCGG") == pytest.approx(1.0)
    # No G, so the expected count is zero and the ratio is undefined.
    assert cpg_observed_expected("AAAA") is None
    assert gc_skew("GGGG") == pytest.approx(1.0)
    assert gc_skew("CCCC") == pytest.approx(-1.0)
    assert gc_skew("GGCC") == pytest.approx(0.0)
    assert gc_skew("AAAA") is None
    # S is G or C, but the skew needs to know which. It stays out.
    assert gc_skew("SSSS") is None


def test_windows_are_one_based_and_do_not_run_past_the_end(tmp_path: Path) -> None:
    fasta = tmp_path / "short.fasta"
    fasta.write_text(">piece\nAAACCCGGGT\n", encoding="utf-8")

    rows = gc_windows(fasta, window=4, step=4)

    assert [(row.start, row.end, round(row.gc_percent, 2), row.cpg_oe) for row in rows] == [
        (1, 4, 25.0, None),
        (5, 8, 100.0, 1.0),
    ]
    # The last two bases are not a full window of 4.
    assert all(row.end - row.start + 1 == 4 for row in rows)


def test_skew_cumulative_restarts_per_record(tmp_path: Path) -> None:
    fasta = tmp_path / "two.fasta"
    fasta.write_text(">one\nGGGGCCCC\n>two\nGGGG\n", encoding="utf-8")

    rows = skew_windows(fasta, window=4, step=4)

    assert [(row.id, row.gc_skew, row.cumulative_skew) for row in rows] == [
        ("one", pytest.approx(1.0), pytest.approx(1.0)),
        ("one", pytest.approx(-1.0), pytest.approx(0.0)),
        ("two", pytest.approx(1.0), pytest.approx(1.0)),
    ]


def test_window_and_step_must_be_positive(tmp_path: Path) -> None:
    fasta = tmp_path / "one.fasta"
    fasta.write_text(">a\nATGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="window"):
        gc_windows(fasta, window=0, step=1)
    with pytest.raises(ValueError, match="step"):
        skew_windows(fasta, window=2, step=0)


def test_sequence_shorter_than_the_window_adds_no_rows(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">a\nATGC\n", encoding="utf-8")
    assert gc_windows(fasta, window=10, step=1) == []
    assert skew_windows(fasta, window=10, step=1) == []
