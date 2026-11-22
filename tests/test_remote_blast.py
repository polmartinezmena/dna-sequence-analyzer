import json

import pytest

from dna_sequence_analyzer import remote_blast as blast


def test_submit_wait_and_retrieve(monkeypatch):
    clock = [100]
    calls = []
    replies = iter(["RID = MOCK123\nRTOE = 4", "Status=WAITING", json.dumps({
        "BlastOutput2": [{"report": {"results": {"search": {"query_id": "test", "hits": [{
            "description": [{"accession": "test-subject"}],
            "hsps": [{"evalue": 0.01, "bit_score": 10, "identity": 8, "align_len": 10}]
        }]}}}}]
    })])

    def request(**params):
        calls.append(params)
        return next(replies)

    monkeypatch.setattr(blast, "_request", request)
    monkeypatch.setattr(blast, "monotonic", lambda: clock[0])
    job = blast.submit("test-query")
    with pytest.raises(ValueError):
        blast.retrieve(job)
    clock[0] += 60
    assert blast.retrieve(job) is None
    clock[0] += 60
    assert blast.retrieve(job)[0].identity_percent == 80
    assert [call["CMD"] for call in calls] == ["Put", "Get", "Get"]


def test_requires_configured_contact(monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(ValueError, match="NCBI_EMAIL"):
        blast.submit("test-query")
