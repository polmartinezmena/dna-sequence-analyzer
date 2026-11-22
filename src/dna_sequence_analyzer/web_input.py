"""Turn text from the web form into the records the rest of the package scores.

The Streamlit app should not reimplement GC, windows or skew. It only
has to decide what the visitor pasted, then hand real SeqRecord objects
to the functions that already exist.
"""

import csv
import io
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dna_sequence_analyzer.gc import SequenceStats, analyze_records

_BASES = frozenset("ACGTU")


def parse_sequence_text(text: str) -> list[SeqRecord]:
    """Read pasted FASTA, or one raw sequence when there is no header.

    A blank line before the first ``>`` is ignored. Whitespace inside a
    raw sequence is ignored, so a sequence copied with spaces still
    works. FASTA is parsed by the same rules as a file: each record
    starts with ``>``.

    Raises:
        ValueError: The text is empty, the FASTA has no records, or a
            raw paste contains no A, C, G, T or U.
    """
    cleaned = text.strip()
    if not cleaned:
        raise ValueError("Pega un FASTA o una secuencia de bases.")
    if looks_like_fasta(cleaned):
        from Bio import SeqIO

        records = list(SeqIO.parse(io.StringIO(cleaned), "fasta"))
        if not records:
            raise ValueError(
                "No encuentro registros FASTA. Cada secuencia tiene que empezar por '>'."
            )
        return records

    sequence = "".join(cleaned.split())
    if not any(base in _BASES for base in sequence.upper()):
        raise ValueError(
            "No veo bases A, C, G, T o U. Revisa que hayas pegado una secuencia y no un texto cualquiera."
        )
    return [SeqRecord(Seq(sequence), id="secuencia", description="secuencia")]


def looks_like_fasta(text: str) -> bool:
    """True when the first non-empty line starts with ``>``."""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped.startswith(">")
    return False


def text_from_upload(data: bytes) -> str:
    """Decode an uploaded file. A FASTA file is plain text."""
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError(
            "No he podido leer el archivo como texto. Un FASTA es un archivo de texto, no un binario."
        ) from exc


def describe_records(records: list[SeqRecord]) -> list[SequenceStats]:
    """GC table for records that came from paste or upload."""
    return analyze_records(records)


def stats_csv(stats: list[SequenceStats]) -> str:
    """CSV with the same columns as ``dna-analyzer gc``."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(["id", "length", "gc_percent", "ambiguous", "type", "flag"])
    for item in stats:
        gc_text = "" if item.gc_percent is None else f"{item.gc_percent:.2f}"
        writer.writerow(
            [
                item.id,
                item.length,
                gc_text,
                item.ambiguous_bases,
                item.nucleic_acid,
                item.gc_flag or "",
            ]
        )
    return buffer.getvalue()


def record_labels(records: list[SeqRecord]) -> list[str]:
    """Labels for a sequence selector. Repeated ids get a counter."""
    totals: dict[str, int] = {}
    for record in records:
        totals[record.id] = totals.get(record.id, 0) + 1
    seen: dict[str, int] = {}
    labels: list[str] = []
    for record in records:
        seen[record.id] = seen.get(record.id, 0) + 1
        if totals[record.id] == 1:
            labels.append(record.id)
        else:
            labels.append(f"{record.id} ({seen[record.id]})")
    return labels


def bundled_example_fasta() -> Path:
    """Find ``data/examples.fasta`` from the repository checkout.

    Streamlit Cloud and ``dna-analyzer web`` both run with the repo
    files available. An installed copy that dropped the ``data`` folder
    gets a clear error instead of a missing-file traceback.
    """
    import sysconfig

    package_dir = Path(__file__).resolve().parent
    candidates = [
        Path(sysconfig.get_path("data")) / "share" / "dna-sequence-analyzer" / "examples.fasta",
        package_dir.parents[1] / "data" / "examples.fasta",
        package_dir.parents[2] / "data" / "examples.fasta",
        Path.cwd() / "data" / "examples.fasta",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        "No encuentro data/examples.fasta. Abre la app desde la carpeta del repositorio."
    )


def explain_tm_error(exc: ValueError) -> str:
    """Spanish message for the Tm form. The library messages are English."""
    message = str(exc)
    if "at least 2" in message or "Sequence is empty" in message:
        return "Hacen falta al menos 2 bases para calcular la Tm."
    if "not removed" in message or "unambiguous" in message:
        return (
            "La Tm solo acepta A, C, G, T o U. No quito las letras ambiguas, "
            "porque eso calcularía otra secuencia."
        )
    if "cannot be negative" in message:
        return "Las concentraciones no pueden ser negativas."
    if "dnac1" in message:
        return (
            "La concentración del oligo tiene que dejar una cantidad positiva "
            "de hebra en la fórmula. Prueba con un valor mayor que 0."
        )
    return "No he podido calcular la Tm con estos datos. Revisa la secuencia y las concentraciones."
