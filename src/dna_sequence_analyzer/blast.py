"""Sequential NCBI BLAST URL API client with a local result cache."""

import json
import math
import os
import re
import sqlite3
from contextlib import closing, contextmanager
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from tempfile import NamedTemporaryFile
from threading import Lock
from time import monotonic as _clock, sleep as _sleep, time as _wall_time
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

_ENDPOINT = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_TOOL = "dna-sequence-analyzer"
_GATE = Lock()
_LAST_CONTACT = None
_DATABASES = {"blastn": ("nt", "refseq_rna"), "blastp": ("nr",)}


@dataclass(frozen=True)
class BlastOptions:
    program: str = "blastn"
    database: str = "nt"
    evalue: float = 0.01
    max_hits: int = 20

    def __post_init__(self):
        if self.database not in _DATABASES.get(self.program, ()):
            raise ValueError("El programa y la base de datos BLAST no son compatibles.")
        if not math.isfinite(self.evalue) or self.evalue <= 0 or self.max_hits < 1:
            raise ValueError("El valor E y el máximo de resultados deben ser positivos.")


@dataclass(frozen=True)
class BlastHit:
    accession: str
    description: str
    identity: float
    coverage: float
    evalue: float
    bitscore: float


def _request(params, timeout):
    encoded = urlencode(params)
    if params["CMD"] == "Put":
        request = Request(_ENDPOINT, data=encoded.encode("ascii"))
    else:
        request = Request(f"{_ENDPOINT}?{encoded}")
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except URLError as error:
        raise OSError("NCBI BLAST no respondió a la consulta.") from error


def _field(response, name):
    field = re.search(rf"^\s*{name}\s*=\s*([A-Za-z0-9_-]+)\s*$", response, re.MULTILINE)
    return field.group(1) if field else None


def _covered_bases(hsps):
    end = total = 0
    intervals = sorted(sorted((hsp["query_from"], hsp["query_to"])) for hsp in hsps)
    for start, stop in intervals:
        total += max(0, stop - max(start - 1, end))
        end = max(end, stop)
    return total


def _parse_report(response, max_hits):
    try:
        reports = json.loads(response)["BlastOutput2"]
        query = reports[0]["report"]["results"]["search"]
        hits = []
        for hit in query["hits"][:max_hits]:
            desc = hit["description"][0]
            best = max(hit["hsps"], key=lambda hsp: (hsp["bit_score"], -hsp["evalue"]))
            hits.append(BlastHit(
                desc["accession"], desc.get("title", ""),
                100 * best["identity"] / best["align_len"],
                100 * _covered_bases(hit["hsps"]) / query["query_len"],
                best["evalue"], best["bit_score"],
            ))
        return hits
    except (ValueError, KeyError, TypeError, IndexError, ZeroDivisionError) as error:
        raise ValueError("NCBI devolvió un informe BLAST que no se pudo interpretar.") from error


def _cache_path(query, options, cache_dir):
    default = Path(os.environ.get("LOCALAPPDATA", Path.home() / ".cache"))
    directory = Path(cache_dir) if cache_dir is not None else default / _TOOL / "blast"
    payload = json.dumps(["v1", query, asdict(options)], sort_keys=True)
    return directory / (sha256(payload.encode()).hexdigest() + ".json")


def _store(path, hits):
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as stream:
        temp = Path(stream.name)
        try:
            json.dump([asdict(hit) for hit in hits], stream)
        except BaseException:
            stream.close()
            temp.unlink(missing_ok=True)
            raise
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


@contextmanager
def _request_queue(directory, timeout):
    """Serialize requests across processes sharing a cache directory."""
    directory.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(directory / "requests.sqlite3", timeout=timeout)) as queue:
        try:
            queue.execute("BEGIN IMMEDIATE")
            queue.execute("CREATE TABLE IF NOT EXISTS contact (stamp REAL NOT NULL)")
        except sqlite3.OperationalError as error:
            raise TimeoutError("Otra búsqueda BLAST todavía está en curso.") from error
        try:
            yield queue
        finally:
            queue.commit()


def search(sequence, options, *, email, timeout=300, cache_dir=None, progress=None):
    """Search one sequence, polling a single RID at intervals of at least one minute."""
    query = "".join(sequence.split()).upper()
    if not query or not query.isascii() or not query.isalpha():
        raise ValueError("Introduce una única secuencia con letras de nucleótidos o aminoácidos.")
    if not email or not re.fullmatch(r"[^\s@]+@[^\s@]+", email):
        raise ValueError("Configura un correo válido en NCBI_EMAIL.")
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("El tiempo de espera debe ser positivo.")
    path = _cache_path(query, options, cache_dir)
    deadline = _clock() + timeout
    if not _GATE.acquire(timeout=timeout):
        raise TimeoutError("Otra búsqueda BLAST todavía está en curso.")
    try:
        remaining = deadline - _clock()
        if remaining <= 0:
            raise TimeoutError("Se agotó el tiempo de espera de BLAST.")
        with _request_queue(path.parent, remaining) as queue:
            return _search_locked(query, options, email, path, deadline, queue, progress)
    finally:
        _GATE.release()


def _search_locked(query, options, email, path, deadline, queue, progress):
    try:
        hits = [BlastHit(**row) for row in json.loads(path.read_text(encoding="utf-8"))]
    except (OSError, ValueError, TypeError):
        hits = None
    if hits is not None:
        if progress:
            progress("Resultado recuperado de la caché.")
        return hits

    contacts = {"EMAIL": email, "TOOL": _TOOL}
    row = queue.execute("SELECT stamp FROM contact").fetchone()
    previous_contact = _clock() + min(0, row[0] - _wall_time()) if row else None

    def send(params, delay=0):
        global _LAST_CONTACT
        nonlocal previous_contact
        for stamp in (previous_contact, _LAST_CONTACT):
            if stamp is not None:
                delay = max(delay, stamp + 10 - _clock())
        if _clock() + delay >= deadline:
            raise TimeoutError("Se agotó el tiempo de espera de BLAST.")
        if delay:
            _sleep(delay)
        _LAST_CONTACT = previous_contact = _clock()
        queue.execute("DELETE FROM contact")
        queue.execute("INSERT INTO contact VALUES (?)", (_wall_time(),))
        return _request({**params, **contacts}, deadline - _clock())

    if progress:
        progress("Enviando la consulta a NCBI.")
    response = send({
        "CMD": "Put", "QUERY": query, "PROGRAM": options.program, "DATABASE": options.database,
        "EXPECT": options.evalue, "HITLIST_SIZE": options.max_hits,
    })
    rid = _field(response, "RID")
    if not rid:
        raise ValueError("NCBI no devolvió un identificador de búsqueda.")
    if progress:
        progress(f"NCBI está procesando la búsqueda ({rid}).")
    delay = max(60, int(_field(response, "RTOE") or 0))
    while True:
        status_info = send({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}, delay)
        status = _field(status_info, "Status")
        if status != "WAITING":
            break
        delay = 60
    if status != "READY":
        raise ValueError("NCBI no pudo completar la búsqueda BLAST.")
    hits = []
    if _field(status_info, "ThereAreHits") != "no":
        if progress:
            progress("Descargando las coincidencias.")
        report = send({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2_S"}, 10)
        hits = _parse_report(report, options.max_hits)
    _store(path, hits)
    return hits
