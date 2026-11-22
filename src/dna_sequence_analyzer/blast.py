"""Summarize alignments from saved BLAST reports."""

from dataclasses import dataclass
from typing import TextIO

from Bio import SearchIO


@dataclass(frozen=True)
class BlastHit:
    query: str
    subject: str
    aligned_columns: int
    identity_percent: float
    evalue: float
    bit_score: float


def summarize_blast(handle: TextIO, report_format: str = "blast-xml") -> list[BlastHit]:
    """Report the lowest E-value segment per query and subject."""
    rows = []
    for query in SearchIO.parse(handle, report_format):
        for hit in query:
            best = min(hit.hsps, key=lambda hsp: (hsp.evalue, -hsp.bitscore))
            identity = (
                best.ident_pct if report_format == "blast-tab"
                else 100 * best.ident_num / best.aln_span
            )
            rows.append(BlastHit(
                query.id, hit.id, best.aln_span, identity, best.evalue, best.bitscore
            ))
    return rows
