"""Tests for the ``dna-analyzer gc`` command."""

from pathlib import Path

from typer.testing import CliRunner

from dna_sequence_analyzer.cli import app

runner = CliRunner()


def test_gc_command_prints_an_aligned_table(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">ab example\nATGC\n>c\nGGG\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout == (
        "id  length  gc_percent\n"
        "--  ------  ----------\n"
        "ab       4       50.00\n"
        "c        3      100.00\n"
    )


def test_gc_command_can_print_csv(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">ab\nATGC\n>weird,id\nAAAN\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta), "--csv"])

    assert result.exit_code == 0
    # A comma inside an id is quoted. The undefined-looking N is excluded,
    # so AAA N is 0% GC, not an empty value. Empty CSV cells are reserved
    # for a percentage that cannot be calculated.
    assert result.stdout == ("id,length,gc_percent\n" "ab,4,50.00\n" '"weird,id",4,0.00\n')


def test_gc_command_shows_na_when_gc_is_undefined(tmp_path: Path) -> None:
    fasta = tmp_path / "unknown.fasta"
    fasta.write_text(">only_n\nNNNN\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert "NA" in result.stdout
    assert "only_n" in result.stdout
    assert "no bases that can be counted" in result.stderr


def test_csv_leaves_undefined_gc_blank(tmp_path: Path) -> None:
    fasta = tmp_path / "unknown.fasta"
    fasta.write_text(">only_n\nNNNN\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta), "--csv"])

    assert result.exit_code == 0
    assert result.stdout == "id,length,gc_percent\nonly_n,4,\n"


def test_directory_is_a_clear_error(tmp_path: Path) -> None:
    result = runner.invoke(app, ["gc", str(tmp_path)])

    assert result.exit_code == 1
    assert "not a file" in result.stderr


def test_missing_file_is_a_clear_error(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fasta"

    result = runner.invoke(app, ["gc", str(missing)])

    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert "not found" in result.stderr


def test_empty_file_is_a_clear_error(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 1
    assert "empty" in result.stderr


def test_invalid_fasta_is_a_clear_error(tmp_path: Path) -> None:
    fasta = tmp_path / "notes.txt"
    fasta.write_text("this is not a fasta file\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 1
    assert "FASTA" in result.stderr


def test_help_lists_the_gc_command() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "gc" in result.stdout
    assert "GC" in result.stdout


def test_version_is_printed() -> None:
    result = runner.invoke(app, ["--version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "dna-analyzer 0.1.0"


def test_example_file_command(tmp_path: Path) -> None:
    example = Path(__file__).resolve().parents[1] / "data" / "examples.fasta"

    result = runner.invoke(app, ["gc", str(example)])

    assert result.exit_code == 0
    assert "NM_000207.3" in result.stdout
    assert "L29345.1" in result.stdout
    assert "465" in result.stdout
    assert "922" in result.stdout
    assert "63.87" in result.stdout
    assert "36.23" in result.stdout
