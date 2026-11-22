from pathlib import Path

import pytest

from dna_sequence_analyzer import blast


FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(("program", "database", "megablast", "expected"), [
    ("blastn", "nt", True, "on"), ("blastn", "nt", False, None), ("blastp", "nr", True, None),
])
def test_submission_polling_results_and_cache(monkeypatch, tmp_path, program, database, megablast, expected):
    monkeypatch.setattr(blast, "_LAST_CONTACT", None)
    now = [0.0]
    requests = []
    files = iter(("blast_put.txt", "blast_waiting.txt", "blast_ready.txt", "blast_result.json"))

    def fake_request(params, timeout):
        requests.append((now[0], params.copy()))
        assert timeout > 0
        return (FIXTURES / next(files)).read_text(encoding="utf-8")

    monkeypatch.setattr(blast, "_request", fake_request)
    monkeypatch.setattr(blast, "_clock", lambda: now[0])
    monkeypatch.setattr(blast, "_sleep", lambda seconds: now.__setitem__(0, now[0] + seconds))
    options = blast.BlastOptions(program=program, database=database, max_hits=5, megablast=megablast)
    query = "MKAA" if program == "blastp" else "ACGT"
    hits = blast.search(query, options, email="test@example.org", cache_dir=tmp_path)
    assert [(h.accession, h.identity, h.coverage, h.evalue, h.bitscore) for h in hits] == [
        ("TEST.1", 90.0, 80.0, 1e-20, 200.0),
    ]
    assert [stamp for stamp, _ in requests] == [0.0, 60.0, 120.0, 130.0]
    assert [params["CMD"] for _, params in requests] == ["Put", "Get", "Get", "Get"]
    assert requests[0][1]["PROGRAM"] == program
    assert requests[0][1].get("MEGABLAST") == expected
    assert all(params["EMAIL"] and params["TOOL"] for _, params in requests)
    assert blast.search(query, options, email="test@example.org", cache_dir=tmp_path) == hits
    assert len(requests) == 4
    assert len(list(tmp_path.glob("*.json"))) == 1


def test_invalid_program_database_combinations():
    with pytest.raises(ValueError, match="base de datos"):
        blast.BlastOptions(program="blastp", database="nt")


def test_timeout_stops_polling(monkeypatch, tmp_path):
    monkeypatch.setattr(blast, "_LAST_CONTACT", None)
    monkeypatch.setattr(blast, "_clock", lambda: 0.0)
    monkeypatch.setattr(blast, "_request", lambda *args: (FIXTURES / "blast_put.txt").read_text())
    with pytest.raises(TimeoutError):
        blast.search("ATGC", blast.BlastOptions(), email="test@example.org",
                     cache_dir=tmp_path, timeout=30)
    assert not list(tmp_path.glob("*.json"))


def test_cli_outputs_only_csv(monkeypatch, tmp_path):
    from typer.testing import CliRunner

    from dna_sequence_analyzer import cli

    fasta = tmp_path / "query.fasta"
    fasta.write_text(";comentario\n>demo\nMKAA\n", encoding="ascii")
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")

    def fake_search(seq, options, **metadata):
        assert seq == "MKAA" and options.database == "nr"
        assert metadata["email"] == "test@example.org"
        return [blast.BlastHit("TEST.1", "Control", 90, 80, 1e-20, 200)]

    monkeypatch.setattr(cli, "search", fake_search)
    result = CliRunner().invoke(cli.app, ["blast", str(fasta), "--program", "blastp"])
    assert result.exit_code == 0, result.output
    assert "TEST.1,Control,90,80,1e-20,200" in result.stdout


def test_megablast_cache_key(tmp_path):
    from dataclasses import replace

    options = blast.BlastOptions()
    assert blast._cache_path("ACGT", options, tmp_path) != blast._cache_path(
        "ACGT", replace(options, megablast=False), tmp_path
    )


@pytest.mark.parametrize(("last", "now", "expected"), [(None, 10, 0), (10, 11.1, 299), (10, 310, 0)])
def test_blast_session_wait(last, now, expected):
    assert blast.session_wait_seconds(last, now) == expected


@pytest.mark.parametrize(("flags", "enabled"), [([], True), (["--no-megablast"], False), (["--megablast"], True)])
def test_cli_passes_megablast_option(monkeypatch, tmp_path, flags, enabled):
    from typer.testing import CliRunner

    from dna_sequence_analyzer import cli

    fasta = tmp_path / "query.fa"
    fasta.write_text(">demo\nACGT\n", encoding="ascii")

    def fake_search(seq, options, **kwargs):
        assert options.megablast is enabled
        return []

    monkeypatch.setattr(cli, "search", fake_search)
    result = CliRunner().invoke(cli.app, ["blast", str(fasta), *flags])
    assert result.exit_code == 0, result.output
