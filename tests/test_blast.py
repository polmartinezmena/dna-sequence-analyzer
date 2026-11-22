from io import StringIO
from typer.testing import CliRunner

from dna_sequence_analyzer.blast import summarize_blast
from dna_sequence_analyzer.cli import app


def test_blast_table_keeps_best_segment_per_subject_and_query():
    report = StringIO(
        "query1\tsubject1\t90\t10\t1\t0\t1\t10\t11\t20\t1e-3\t20\n"
        "query1\tsubject1\t80\t20\t4\t0\t1\t20\t41\t22\t1e-6\t40\n"
        "query2\tsubject1\t100\t10\t0\t0\t1\t10\t1\t10\t1e-5\t35\n"
    )
    hits = summarize_blast(report, "blast-tab")
    assert [(hit.query, hit.subject, hit.identity_percent, hit.evalue) for hit in hits] == [
        ("query1", "subject1", 80.0, 1e-6),
        ("query2", "subject1", 100.0, 1e-5),
    ]


def test_empty_tabular_report_has_no_hits():
    assert summarize_blast(StringIO(""), "blast-tab") == []


def test_cli_exports_blast_table_as_csv(tmp_path):
    report = tmp_path / "report.tsv"
    report.write_text("q\ts\t90\t10\t1\t0\t1\t10\t11\t20\t1e-3\t20\n", encoding="ascii")
    result = CliRunner().invoke(app, ["blast-report", str(report), "--tabular"])
    assert result.exit_code == 0, result.output
    assert "q,s,10,90.00,0.001,20.0" in result.stdout
