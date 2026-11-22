"""Tests for the ``dna-analyzer gc`` command."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from dna_sequence_analyzer.cli import app

runner = CliRunner()


@pytest.mark.parametrize("command", ["gc", "gc-window", "skew", "translate", "tm", "orfs", "blast", "web"])
def test_command_help_exits_successfully(command):
    result = runner.invoke(app, [command, "--help"], prog_name="dna-analyzer")
    assert result.exit_code == 0
    assert f"dna-analyzer {command}" in result.output


def test_gc_command_prints_an_aligned_table(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">ab example\nATGC\n>c\nGGG\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout == (
        "id  length  gc_percent  ambiguous         type     flag\n"
        "--  ------  ----------  ---------  -----------  -------\n"
        "ab       4       50.00          0          dna\n"
        "c        3      100.00          0  unspecified  high-gc\n"
    )
    assert "high-gc" in result.stderr


def test_gc_command_can_print_csv(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">ab\nATGC\n>weird,id\nAAAN\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta), "--csv"])

    assert result.exit_code == 0
    # A comma inside an id is quoted. The undefined-looking N is excluded,
    # so AAA N is 0% GC, not an empty value. Empty CSV cells are reserved
    # for a percentage that cannot be calculated.
    assert result.stdout == (
        "id,length,gc_percent,ambiguous,type,flag\n"
        "ab,4,50.00,0,dna,\n"
        '"weird,id",4,0.00,1,unspecified,low-gc\n'
    )


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
    assert result.stdout == "id,length,gc_percent,ambiguous,type,flag\nonly_n,4,,4,unspecified,\n"


def test_directory_is_a_clear_error(tmp_path: Path) -> None:
    result = runner.invoke(app, ["gc", str(tmp_path)])

    assert result.exit_code == 1
    assert "no es un archivo" in result.stderr


def test_missing_file_is_a_clear_error(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fasta"

    result = runner.invoke(app, ["gc", str(missing)])

    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert "No se encontró" in result.stderr


def test_empty_file_is_a_clear_error(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 1
    assert "está vacío" in result.stderr


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
    assert result.stdout.strip() == "dna-analyzer 1.0.0"


def test_tm_command_prints_nearest_neighbor_temperature(tmp_path: Path) -> None:
    from dna_sequence_analyzer.tm import melting_temperature

    fasta = tmp_path / "oligos.fasta"
    fasta.write_text(">ab\nATATATAT\n>block\nAAAATTTT\n", encoding="utf-8")

    result = runner.invoke(app, ["tm", str(fasta)])

    alternating = f"{melting_temperature('ATATATAT'):.2f}"
    blocked = f"{melting_temperature('AAAATTTT'):.2f}"
    assert result.exit_code == 0
    assert result.stdout.splitlines()[0] == "id     length  tm_celsius"
    assert alternating in result.stdout
    assert blocked in result.stdout
    assert alternating != blocked
    assert "short DNA oligos" not in result.stderr


def test_tm_sequence_option_and_salt() -> None:
    from dna_sequence_analyzer.tm import MeltingConditions, melting_temperature

    plain = runner.invoke(app, ["tm", "--sequence", "ATGCATGC"])
    salty = runner.invoke(app, ["tm", "--sequence", "ATGCATGC", "--na", "200"])

    assert plain.exit_code == 0
    assert f"{melting_temperature('ATGCATGC'):.2f}" in plain.stdout
    assert salty.exit_code == 0
    changed = melting_temperature("ATGCATGC", MeltingConditions(na_mM=200))
    assert f"{changed:.2f}" in salty.stdout
    assert plain.stdout != salty.stdout


def test_tm_requires_a_file_or_a_sequence(tmp_path: Path) -> None:
    missing = runner.invoke(app, ["tm"])
    both = runner.invoke(app, ["tm", str(tmp_path / "x.fasta"), "--sequence", "ATGC"])

    assert missing.exit_code == 1
    assert "archivo FASTA o --sequence" in missing.stderr
    assert both.exit_code == 1
    assert "única entrada" in both.stderr


def test_tm_missing_file_mentions_sequence_option(tmp_path: Path) -> None:
    result = runner.invoke(app, ["tm", str(tmp_path / "missing.fasta")])

    assert result.exit_code == 1
    assert "No se encontró" in result.stderr
    assert "--sequence" in result.stderr


def test_tm_csv_and_undefined_records(tmp_path: Path) -> None:
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(">ok\nATGC\n>n\nAAAN\n", encoding="utf-8")

    result = runner.invoke(app, ["tm", str(fasta), "--csv"])

    assert result.exit_code == 0
    assert result.stdout.splitlines()[0] == "id,length,tm_celsius"
    assert result.stdout.splitlines()[2] == "n,4,"
    assert "Ambiguous bases are not removed" in result.stderr


def test_tm_warns_that_long_records_are_not_primers(tmp_path: Path) -> None:
    fasta = tmp_path / "long.fasta"
    fasta.write_text(">transcript\n" + ("ATGC" * 20) + "\n", encoding="utf-8")

    result = runner.invoke(app, ["tm", str(fasta)])

    assert result.exit_code == 0
    assert "not for long molecules" in result.stderr
    assert "transcript" in result.stdout


def test_tm_rejects_a_negative_salt() -> None:
    result = runner.invoke(app, ["tm", "--sequence", "ATGC", "--na", "-5"])

    assert result.exit_code == 1
    assert "cannot be negative" in result.stderr


def test_gc_flags_rna_and_reports_ambiguous_bases(tmp_path: Path) -> None:
    fasta = tmp_path / "rna.fasta"
    fasta.write_text(">oligo\naugcnn\n", encoding="utf-8")

    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert "rna" in result.stdout
    assert "2" in result.stdout
    assert "G-U wobble" in result.stderr


def test_gc_window_and_skew_commands(tmp_path: Path) -> None:
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(">one\nGGGGCCCC\n", encoding="utf-8")

    gc_result = runner.invoke(app, ["gc-window", str(fasta), "--window", "4", "--step", "4"])
    skew_result = runner.invoke(app, ["skew", str(fasta), "--window", "4", "--step", "4", "--csv"])

    assert gc_result.exit_code == 0
    assert "cpg_oe" in gc_result.stdout
    assert "100.00" in gc_result.stdout
    assert skew_result.exit_code == 0
    assert skew_result.stdout.splitlines()[0] == "id,start,end,gc_skew,cumulative_skew"
    assert "1.0000" in skew_result.stdout
    assert "-1.0000" in skew_result.stdout


def test_gc_window_explains_when_nothing_fits(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">a\nATGC\n", encoding="utf-8")

    result = runner.invoke(app, ["gc-window", str(fasta), "--window", "10", "--step", "1"])

    assert result.exit_code == 0
    assert "shorter than --window" in result.stderr


def test_tm_warns_when_uracil_is_read_as_thymine() -> None:
    result = runner.invoke(app, ["tm", "--sequence", "AUGCAUGC"])

    assert result.exit_code == 0
    assert "G-U wobble" in result.stderr


def test_plot_without_matplotlib_is_a_clear_error(tmp_path: Path, monkeypatch) -> None:
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(">one\nGGGGCCCC\n", encoding="utf-8")
    import builtins

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ImportError("matplotlib blocked for this test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    result = runner.invoke(
        app,
        ["skew", str(fasta), "--window", "4", "--step", "4", "--plot", str(tmp_path / "skew.png")],
    )

    assert result.exit_code == 1
    assert "matplotlib" in result.stderr


def test_skew_plot_writes_a_png(tmp_path: Path) -> None:
    pytest = __import__("pytest")
    pytest.importorskip("matplotlib")
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(">one\nGGGGCCCC\n>two\nCCCCCCCC\n", encoding="utf-8")
    png = tmp_path / "skew.png"

    result = runner.invoke(
        app,
        ["skew", str(fasta), "--window", "4", "--step", "4", "--plot", str(png)],
    )

    assert result.exit_code == 0
    assert png.is_file()
    assert png.stat().st_size > 0
    assert f"Archivo guardado en {png}" in result.stderr


def test_help_lists_the_tm_command() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "tm" in result.stdout


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
