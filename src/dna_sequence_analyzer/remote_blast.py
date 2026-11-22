"""Interactive NCBI BLAST URL API client for a single application process."""

from dataclasses import dataclass
import json
import os
import re
from threading import Lock
from time import monotonic
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .blast import BlastHit

_lock = Lock()
_last_request = float("-inf")


@dataclass
class SearchJob:
    rid: str
    next_poll: float


def _request(**params):
    global _last_request
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        raise ValueError("Configura NCBI_EMAIL en el entorno que ejecuta Streamlit.")
    with _lock:
        if monotonic() - _last_request < 10:
            raise ValueError("Espera diez segundos antes de contactar de nuevo con NCBI.")
        _last_request = monotonic()
    body = urlencode(dict(params, email=email, tool="dna-sequence-analyzer")).encode()
    request = Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=body)
    with urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def submit(query: str) -> SearchJob:
    """Submit one accession or nucleotide sequence, returning the NCBI RID."""
    reply = _request(CMD="Put", PROGRAM="blastn", DATABASE="core_nt", QUERY=query, HITLIST_SIZE=20)
    match = re.search(r"^\s*RID\s*=\s*([A-Za-z0-9-]+)\s*$", reply, re.M)
    if not match:
        raise ValueError("NCBI no devolvió un identificador de búsqueda.")
    estimate = re.search(r"RTOE\s*=\s*(\d+)", reply)
    return SearchJob(match[1], monotonic() + max(60, int(estimate[1]) if estimate else 60))


def retrieve(job: SearchJob) -> list[BlastHit] | None:
    """Return None while queued and an empty list for a completed search without hits."""
    if monotonic() < job.next_poll:
        raise ValueError("La búsqueda está en cola. Espera al siguiente intervalo de consulta.")
    job.next_poll = monotonic() + 60
    reply = _request(CMD="Get", RID=job.rid, FORMAT_TYPE="JSON2_S")
    if re.search(r"Status\s*=\s*WAITING", reply):
        return None
    try:
        payload = json.loads(reply)["BlastOutput2"]
        entry = payload[0] if isinstance(payload, list) else payload
        search = entry["report"]["results"]["search"]
        hits = []
        for hit in search["hits"]:
            best = min(hit["hsps"], key=lambda h: (h["evalue"], -h["bit_score"]))
            hits.append(BlastHit(search["query_id"], hit["description"][0]["accession"],
                                 best["align_len"], 100 * best["identity"] / best["align_len"],
                                 best["evalue"], best["bit_score"]))
        return hits
    except (KeyError, IndexError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"NCBI no devolvió un informe válido para {job.rid}.") from exc
