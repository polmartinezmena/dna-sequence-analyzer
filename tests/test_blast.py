from pathlib import Path

import pytest

from dna_sequence_analyzer import blast


FIXTURES = Path(__file__).parent / "fixtures"


def test_submission_polling_results_and_cache(monkeypatch, tmp_path):
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
    options = blast.BlastOptions(program="blastp", database="nr", max_hits=5)
    hits = blast.search("MKAA", options, email="test@example.org", cache_dir=tmp_path)
    assert [(h.accession, h.identity, h.coverage, h.evalue, h.bitscore) for h in hits] == [
        ("TEST.1", 90.0, 80.0, 1e-20, 200.0),
    ]
    assert [stamp for stamp, _ in requests] == [0.0, 60.0, 120.0, 130.0]
    assert [params["CMD"] for _, params in requests] == ["Put", "Get", "Get", "Get"]
    assert requests[0][1]["PROGRAM"] == "blastp"
    assert all(params["EMAIL"] and params["TOOL"] for _, params in requests)
    assert blast.search("MKAA", options, email="test@example.org", cache_dir=tmp_path) == hits
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
    fasta.write_text(">demo\nMKAA\n", encoding="ascii")
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")

    def fake_search(seq, options, **metadata):
        assert seq == "MKAA" and options.database == "nr"
        assert metadata["email"] == "test@example.org"
        return [blast.BlastHit("TEST.1", "Control", 90, 80, 1e-20, 200)]

    monkeypatch.setattr(cli, "search", fake_search)
    result = CliRunner().invoke(cli.app, ["blast", str(fasta), "--program", "blastp"])
    assert result.exit_code == 0, result.output
    assert "TEST.1,Control,90,80,1e-20,200" in result.stdout
