"""Tests for FASTA reading and the GC-content calculation."""

from pathlib import Path

import pytest

from dna_sequence_analyzer.gc import analyze_fasta, gc_percent

EXAMPLE_FASTA = Path(__file__).resolve().parents[1] / "data" / "examples.fasta"


def write_fasta(path: Path, records: list[tuple[str, str]], *, newline: str = "\n") -> None:
    """Write a tiny FASTA file. ``records`` is a list of (header, sequence)."""
    lines: list[str] = []
    for header, sequence in records:
        lines.append(f">{header}")
        lines.append(sequence)
    path.write_text(newline.join(lines) + newline, encoding="utf-8")


def test_gc_percent_of_a_known_sequence() -> None:
    # 2 GC bases out of 4 unambiguous bases.
    assert gc_percent("ATGC") == pytest.approx(50.0)
    assert gc_percent("GGG") == pytest.approx(100.0)
    # Zero is a real percentage. It is not the same as "undefined".
    assert gc_percent("AAA") == pytest.approx(0.0)


def test_gc_percent_ignores_letter_case() -> None:
    assert gc_percent("atgc") == pytest.approx(50.0)
    assert gc_percent("AtGc") == pytest.approx(50.0)


def test_ambiguous_bases_are_excluded_from_the_percentage() -> None:
    # The four Ns are unknown, so they must not dilute the 50% of ATGC.
    assert gc_percent("ATGCNNNN") == pytest.approx(50.0)
    assert gc_percent("atgcnn") == pytest.approx(50.0)
    # S is "G or C": either base is GC. W is "A or T": neither base is GC.
    assert gc_percent("SSSS") == pytest.approx(100.0)
    assert gc_percent("WWWW") == pytest.approx(0.0)
    assert gc_percent("SW") == pytest.approx(50.0)
    # R is "A or G", which may or may not be GC, so the R is left out.
    assert gc_percent("ATGR") == pytest.approx(100.0 / 3)


def test_undefined_gc_is_none_not_zero() -> None:
    assert gc_percent("") is None
    assert gc_percent("NNNN") is None
    assert gc_percent("----") is None
    assert gc_percent("nry") is None


def test_rna_uracil_counts_as_an_unambiguous_non_gc_base() -> None:
    # U is the RNA form of T. AUGC is half G/C, just like ATGC.
    assert gc_percent("AUGC") == pytest.approx(50.0)
    assert gc_percent("UUU") == pytest.approx(0.0)


def test_whitespace_inside_a_string_is_ignored() -> None:
    assert gc_percent("AT GC") == pytest.approx(50.0)


def test_analyze_fasta_reads_several_records(tmp_path: Path) -> None:
    fasta = tmp_path / "multi.fasta"
    write_fasta(
        fasta,
        [
            ("alpha a first sequence", "ATGC"),
            ("beta", "GGGGAAAA"),
            ("alpha", "CC"),
        ],
    )

    stats = analyze_fasta(fasta)

    assert [item.id for item in stats] == ["alpha", "beta", "alpha"]
    assert [item.name for item in stats] == ["alpha", "beta", "alpha"]
    assert stats[0].description == "alpha a first sequence"
    assert stats[0].length == 4
    assert stats[0].gc_percent == pytest.approx(50.0)
    assert stats[0].gc_bases == 2
    assert stats[0].unambiguous_bases == 4
    assert stats[0].ambiguous_bases == 0
    # 4 G and 4 A: GC% is 50, and the relationship is the documented formula.
    assert stats[1].gc_percent == pytest.approx(100.0 * stats[1].gc_bases / stats[1].unambiguous_bases)
    assert stats[1].gc_bases == 4
    assert stats[2].gc_percent == pytest.approx(100.0)


def test_wrapped_lines_belong_to_one_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "wrapped.fasta"
    fasta.write_text(">wrapped example\nATGC\nATGC\n", encoding="utf-8")

    stats = analyze_fasta(fasta)

    assert len(stats) == 1
    assert stats[0].length == 8
    assert stats[0].gc_percent == pytest.approx(50.0)


def test_lowercase_and_ambiguous_bases_in_a_file(tmp_path: Path) -> None:
    fasta = tmp_path / "soft.fasta"
    write_fasta(fasta, [("masked", "atgcNNNN")])

    stats = analyze_fasta(fasta)

    assert stats[0].length == 8
    assert stats[0].gc_percent == pytest.approx(50.0)
    assert stats[0].unambiguous_bases == 4
    assert stats[0].ambiguous_bases == 4


def test_gaps_count_toward_length_but_not_toward_gc(tmp_path: Path) -> None:
    fasta = tmp_path / "gapped.fasta"
    write_fasta(fasta, [("gapped", "AT-GC")])

    stats = analyze_fasta(fasta)

    assert stats[0].length == 5
    assert stats[0].gc_bases == 2
    assert stats[0].unambiguous_bases == 4
    assert stats[0].ambiguous_bases == 0
    assert stats[0].gc_percent == pytest.approx(50.0)


def test_header_only_record_has_undefined_gc(tmp_path: Path) -> None:
    fasta = tmp_path / "empty_record.fasta"
    fasta.write_text(">empty\n>full\nATGC\n", encoding="utf-8")

    stats = analyze_fasta(fasta)

    assert [(item.id, item.length, item.gc_percent) for item in stats] == [
        ("empty", 0, None),
        ("full", 4, pytest.approx(50.0)),
    ]


def test_all_ambiguous_record(tmp_path: Path) -> None:
    fasta = tmp_path / "unknown.fasta"
    write_fasta(fasta, [("unknown", "NNNN")])

    stats = analyze_fasta(fasta)

    assert stats[0].length == 4
    assert stats[0].gc_percent is None
    assert stats[0].ambiguous_bases == 4
    assert stats[0].unambiguous_bases == 0


def test_missing_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "no-such-file.fasta"
    with pytest.raises(FileNotFoundError, match="No se encontró"):
        analyze_fasta(missing)


def test_directory_path_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no es un archivo"):
        analyze_fasta(tmp_path)


def test_empty_file_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="está vacío"):
        analyze_fasta(fasta)


def test_whitespace_only_file_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "blank.fasta"
    fasta.write_text("\n\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cabecera"):
        analyze_fasta(fasta)


def test_file_without_a_header_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "plain.txt"
    fasta.write_text("ATGC\nATGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="FASTA"):
        analyze_fasta(fasta)


def test_binary_file_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "binary.fasta"
    fasta.write_bytes(b"\xff\xfe not text \x00")
    with pytest.raises(ValueError, match="texto sin formato|FASTA"):
        analyze_fasta(fasta)


def test_windows_newlines(tmp_path: Path) -> None:
    fasta = tmp_path / "crlf.fasta"
    write_fasta(fasta, [("crlf", "ATGC")], newline="\r\n")

    stats = analyze_fasta(fasta)

    assert stats[0].id == "crlf"
    assert stats[0].length == 4
    assert stats[0].gc_percent == pytest.approx(50.0)


def test_example_file_matches_the_published_sequences() -> None:
    """The example file holds two unmodified public NCBI transcripts.

    Insulin mRNA (NM_000207.3) is 465 bp with 297 G/C bases.
    GFP mRNA (L29345.1) is 922 bp with 334 G/C bases.
    Those counts were taken from the records themselves, not from the
    function under test.
    """
    stats = analyze_fasta(EXAMPLE_FASTA)

    assert [(item.id, item.length, item.gc_bases, item.ambiguous_bases) for item in stats] == [
        ("NM_000207.3", 465, 297, 0),
        ("L29345.1", 922, 334, 0),
    ]
    assert stats[0].gc_percent == pytest.approx(100.0 * 297 / 465)
    assert stats[1].gc_percent == pytest.approx(100.0 * 334 / 922)
    assert stats[0].description.startswith("NM_000207.3 Homo sapiens insulin")
    assert "green-fluorescent protein" in stats[1].description


def test_s_and_w_count_as_decided_bases(tmp_path: Path) -> None:
    fasta = tmp_path / "strong_weak.fasta"
    write_fasta(fasta, [("sw", "SSww")])

    stats = analyze_fasta(fasta)

    assert stats[0].length == 4
    assert stats[0].gc_bases == 2
    assert stats[0].unambiguous_bases == 4
    assert stats[0].ambiguous_bases == 0
    assert stats[0].gc_percent == pytest.approx(50.0)


def test_gc_percent_agrees_with_biopython_when_bases_are_countable() -> None:
    """Our percentage is Biopython's fraction, times 100.

    ``ambiguous="remove"`` keeps G, C and S as GC and keeps A, T, U and
    W in the total. Other ambiguity codes are dropped. Biopython
    returns 0 when nothing is left to count; this package returns None
    so the table can show NA instead of a misleading 0%.
    """
    from Bio.SeqUtils import gc_fraction

    sequence = "atgCnnnSww"
    assert gc_percent(sequence) == pytest.approx(100.0 * 3 / 7)
    assert gc_percent(sequence) == pytest.approx(gc_fraction(sequence, ambiguous="remove") * 100)
    assert gc_fraction("NNNN", ambiguous="remove") == 0
    assert gc_percent("NNNN") is None
