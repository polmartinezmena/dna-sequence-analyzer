"""Parse FASTA with an optional comment preamble."""

from io import StringIO
from typing import TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_fasta_text(text: str | TextIO) -> tuple[list[SeqRecord], int]:
    """Skip initial blanks and comments, then parse the records with Biopython."""
    handle = StringIO(text) if isinstance(text, str) else text
    comments = 0
    while True:
        position = handle.tell()
        line = handle.readline()
        if not line:
            return [], comments
        if not line.strip():
            continue
        if line.lstrip().startswith((";", "#")):
            comments += 1
            continue
        if not line.startswith(">"):
            raise ValueError("Cada registro FASTA debe empezar por '>'.")
        handle.seek(position)
        return list(SeqIO.parse(handle, "fasta")), comments
