"""Tests for the text the web app accepts before it calls the calculators."""

from pathlib import Path

import pytest

from dna_sequence_analyzer.cli import streamlit_entrypoint, web_launch_command
from dna_sequence_analyzer.web_input import (
    bundled_example_fasta,
    describe_records,
    explain_tm_error,
    looks_like_fasta,
    parse_sequence_text,
    record_labels,
    stats_csv,
    text_from_upload,
)
from dna_sequence_analyzer.windows import gc_windows, gc_windows_for_records, skew_extrema, skew_windows


def test_raw_paste_becomes_one_sequence() -> None:
    records = parse_sequence_text("at gc\natgc")
    assert len(records) == 1
    assert records[0].id == "secuencia"
    assert str(records[0].seq) == "atgcatgc"
    stats = describe_records(records)
    assert stats[0].length == 8
    assert stats[0].gc_percent == pytest.approx(50.0)
    assert stats[0].nucleic_acid == "dna"


def test_pasted_fasta_keeps_every_record() -> None:
    text = "\n>alfa ejemplo\nATGC\n>beta\nGGGG\n"
    assert looks_like_fasta(text) is True
    records = parse_sequence_text(text)
    assert [record.id for record in records] == ["alfa", "beta"]
    assert describe_records(records)[1].gc_percent == pytest.approx(100.0)


def test_empty_and_non_sequence_text_are_rejected() -> None:
    with pytest.raises(ValueError, match="Pega un FASTA"):
        parse_sequence_text("   \n")
    with pytest.raises(ValueError, match="No veo bases"):
        parse_sequence_text("xyz xyz")
    empty_header = parse_sequence_text(">\n")
    assert len(empty_header) == 1
    assert str(empty_header[0].seq) == ""
    assert describe_records(empty_header)[0].gc_percent is None


def test_upload_must_be_text() -> None:
    assert text_from_upload("ATGC\n".encode("utf-8")) == "ATGC\n"
    with pytest.raises(ValueError, match="texto"):
        text_from_upload(b"\xff\xfe no es utf-8")


def test_repeated_ids_get_distinct_labels() -> None:
    records = parse_sequence_text(">gen\nATGC\n>gen\nGGGG\n>otro\nAAAA\n")
    assert record_labels(records) == ["gen (1)", "gen (2)", "otro"]


def test_csv_matches_the_command_columns() -> None:
    stats = describe_records(parse_sequence_text(">ab\nATGC\n>n\nNNNN\n"))
    assert stats_csv(stats) == (
        "id,length,gc_percent,ambiguous,type,flag\n"
        "ab,4,50.00,0,dna,\n"
        "n,4,,4,unspecified,\n"
    )


def test_example_file_is_the_bundled_fasta() -> None:
    path = bundled_example_fasta()
    assert path.name == "examples.fasta"
    stats = describe_records(parse_sequence_text(path.read_text(encoding="utf-8")))
    assert [item.id for item in stats] == ["NM_000207.3", "L29345.1"]


def test_in_memory_windows_match_the_file_reader(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fasta"
    text = ">one\nGGGGCCCC\n"
    fasta.write_text(text, encoding="utf-8")
    records = parse_sequence_text(text)
    assert gc_windows_for_records(records, 4, 4) == gc_windows(fasta, 4, 4)
    rows = skew_windows(fasta, 4, 4)
    lowest, highest = skew_extrema(rows)
    assert lowest is not None and highest is not None
    assert highest.cumulative_skew == pytest.approx(1.0)
    assert highest.start == 1
    assert lowest.cumulative_skew == pytest.approx(0.0)
    assert lowest.start == 5


def test_tm_errors_are_translated() -> None:
    assert "al menos 2" in explain_tm_error(ValueError("Nearest-neighbor Tm needs at least 2 bases."))
    assert "otra secuencia" in explain_tm_error(ValueError("Ambiguous letters are not removed"))
    assert "negativ" in explain_tm_error(ValueError("na_mM cannot be negative."))


def test_web_command_points_at_the_launcher() -> None:
    entry = streamlit_entrypoint()
    assert entry.name == "streamlit_app.py"
    assert entry.is_file()
    command = web_launch_command(8502)
    assert command[1:4] == ["-m", "streamlit", "run"]
    assert command[-1] == "8502"
    assert command[-2] == "--server.port"


def test_web_command_explains_a_missing_extra(monkeypatch) -> None:
    import importlib.util

    from typer.testing import CliRunner

    from dna_sequence_analyzer.cli import app

    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    result = CliRunner().invoke(app, ["web"])
    assert result.exit_code == 1
    assert "Streamlit" in result.stderr


def test_fasta_preamble_does_not_reach_biopython(tmp_path):
    import warnings

    from dna_sequence_analyzer.gc import read_fasta
    from dna_sequence_analyzer.lab import quality
    from dna_sequence_analyzer.web_input import read_fasta_text

    plain = ">uno\nACGT\n>dos\nGG\n"
    annotated = "\n; comentario\n  # otro\n" + plain
    path = tmp_path / "comentado.fa"
    path.write_text(annotated, encoding="utf-8")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        expected, _ = read_fasta_text(plain)
        records, comments = read_fasta_text(annotated)
        for parsed in (records, read_fasta(path), parse_sequence_text(annotated)):
            assert [(r.id, r.seq) for r in parsed] == [(r.id, r.seq) for r in expected]
    assert ("Aviso", "se ignoraron 2 líneas de comentario antes del primer registro") in quality(
        records, describe_records(records), comments=comments
    )
