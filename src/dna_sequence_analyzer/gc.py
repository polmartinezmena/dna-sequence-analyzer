"""Read FASTA records and calculate their GC content.

GC content is the share of a DNA (or RNA) sequence that is guanine
or cytosine. In a double helix, G pairs with C through three hydrogen
bonds, while A pairs with T through two. A higher GC content therefore
tends to make the duplex harder to melt, which is one reason the number
is used to describe genes, genomes and primers.

The percentage this module returns is::

    GC% = 100 * (G + C + S) / (A + C + G + T + U + S + W)

``S`` means "G or C" and ``W`` means "A or T", so those two letters
still tell us which side of the fraction a base belongs to. Letters
that could fall on either side, such as ``N``, are left out. See
``gc_percent`` for the full rule.
"""

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Bases that must be guanine or cytosine.
# S is the IUPAC "strong" code: the base is G or C, but the call did
# not say which. Either way it forms three hydrogen bonds, so it counts
# as GC.
_GC_BASES = frozenset("GCS")
# Bases that are known and are not guanine or cytosine.
# U is uracil, the RNA equivalent of T. W is the IUPAC "weak" code:
# the base is A or T, so it does not count as GC, but it does count in
# the total.
_AT_BASES = frozenset("ATUW")
# Ambiguity codes that could be either a GC pair or an AT pair.
# They are left out of both sides of the percentage.
#   R = A or G      Y = C or T      K = G or T      M = A or C
#   B = not A       D = not C       H = not G       V = not T
#   N = any base    X = any base (a non-standard letter used like N)
_AMBIGUOUS_BASES = frozenset("RYKMBDHVNX")


@dataclass(frozen=True)
class SequenceStats:
    """Length and GC content of one FASTA record.

    Attributes:
        id: Sequence name. This is the first word of the FASTA header,
            the text that comes right after ``>``.
        description: Full header line, without the leading ``>``.
        length: Number of symbols in the sequence. Ambiguous bases,
            gaps and any other characters are included.
        gc_percent: GC content as a percentage between 0 and 100.
            ``None`` when the record has no unambiguous bases, because
            the percentage would require a division by zero.
        gc_bases: How many symbols count as G or C. That is the letters
            G, C and S.
        unambiguous_bases: How many symbols have a known side in the
            GC fraction: A, C, G, T, U, S and W.
        ambiguous_bases: How many letters were excluded because they
            could be either GC or AT (N, R, Y, K, M, B, D, H, V and X).
    """

    id: str
    description: str
    length: int
    gc_percent: float | None
    gc_bases: int
    unambiguous_bases: int
    ambiguous_bases: int

    @property
    def name(self) -> str:
        """Sequence name. In FASTA this is the same value as ``id``."""
        return self.id


def gc_percent(sequence: str) -> float | None:
    """Return the GC percentage of a DNA or RNA sequence.

    The calculation is::

        GC% = 100 * (G + C + S) / (A + C + G + T + U + S + W)

    Letter case does not matter: ``atgc`` and ``ATGC`` are the same
    sequence. Sequencing files are often written in lowercase, and
    repeat-masking software uses lowercase to mark repetitive DNA.
    The base is still an A, C, G or T, so it still counts.

    Two ambiguous letters are decided enough to keep:

    * ``S`` means G or C ("strong", three hydrogen bonds). It counts
      as GC.
    * ``W`` means A or T ("weak", two hydrogen bonds). It counts in
      the total, and not as GC.

    The other ambiguity codes could be either kind of pair, so they
    are left out of both the numerator and the denominator. An ``N``
    means "base unknown". Counting those positions in the denominator
    would make the GC content look lower than the known sequence
    really is, and counting them as G or C would be a guess. The same
    is true of ``R``, ``Y``, ``K``, ``M``, ``B``, ``D``, ``H``, ``V``
    and of ``X``, which some files use in the same way as ``N``.

    Gap characters (``-``, ``.``), whitespace and any other symbol are
    ignored here as well. They do count toward :attr:`SequenceStats.length`
    when the sequence comes from a FASTA record.

    Returns ``None`` when the sequence contains none of A, C, G, T, U,
    S or W. A sequence of only Ns has an undefined GC content; it is
    not 0%. A sequence of only adenines is defined, and the result is 0%.

    This is the same rule as Biopython's ``Bio.SeqUtils.gc_fraction``
    with ``ambiguous="remove"``. That function returns a fraction
    between 0 and 1, and it returns 0 when nothing can be counted.
    This function returns a percentage between 0 and 100, which is the
    number printed by the command-line tool, and it returns ``None``
    in the undefined case so that "no information" is not printed as
    if the sequence were all A and T.
    """
    gc_bases, at_bases, _ambiguous_bases = _count_bases(sequence)
    unambiguous_bases = gc_bases + at_bases
    if unambiguous_bases == 0:
        return None
    return 100.0 * gc_bases / unambiguous_bases


def analyze_fasta(path: str | Path) -> list[SequenceStats]:
    """Read a FASTA file and describe every record.

    The file may contain one sequence or many. Records are returned in
    the same order as they appear in the file. Wrapped sequence lines
    are joined, which is the normal FASTA layout: a sequence is often
    broken into lines of 60 or 70 characters and those line breaks are
    not part of the molecule.

    Args:
        path: Path to a FASTA file.

    Returns:
        One :class:`SequenceStats` per record. A header with an empty
        sequence is kept: its length is 0 and its GC content is
        ``None``.

    Raises:
        FileNotFoundError: ``path`` does not exist.
        ValueError: ``path`` is not a file, the file is empty, it is
            not text, or it contains no FASTA records. A FASTA record
            has to start with a ``>`` header line.
    """
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not fasta_path.is_file():
        raise ValueError(f"Path is not a file: {fasta_path}")
    if fasta_path.stat().st_size == 0:
        raise ValueError(f"FASTA file is empty: {fasta_path}")

    try:
        with fasta_path.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except UnicodeDecodeError as exc:
        raise ValueError(
            f"Could not read {fasta_path} as text. A FASTA file should be plain text."
        ) from exc
    except ValueError as exc:
        raise ValueError(
            f"Could not parse {fasta_path} as FASTA. "
            "Each record must start with a '>' header line."
        ) from exc

    if not records:
        raise ValueError(
            f"No FASTA records found in {fasta_path}. "
            "Each record must start with a '>' header line followed by sequence."
        )

    return [_stats_from_record(record) for record in records]


def _count_bases(sequence: str) -> tuple[int, int, int]:
    """Count GC symbols, AT symbols and excluded ambiguity codes.

    The first number is G, C and S. The second is A, T, U and W.
    The third is the ambiguity codes that are left out of the
    percentage. The count is case-insensitive. Anything else (gaps,
    spaces, letters outside the nucleic-acid alphabet) is not included.
    """
    gc_bases = 0
    at_bases = 0
    ambiguous_bases = 0
    for base in sequence.upper():
        if base in _GC_BASES:
            gc_bases += 1
        elif base in _AT_BASES:
            at_bases += 1
        elif base in _AMBIGUOUS_BASES:
            ambiguous_bases += 1
    return gc_bases, at_bases, ambiguous_bases


def _stats_from_record(record: SeqRecord) -> SequenceStats:
    """Build a :class:`SequenceStats` from one parsed FASTA record."""
    sequence = str(record.seq)
    gc_bases, at_bases, ambiguous_bases = _count_bases(sequence)
    unambiguous_bases = gc_bases + at_bases
    percent = None if unambiguous_bases == 0 else 100.0 * gc_bases / unambiguous_bases
    return SequenceStats(
        id=record.id,
        description=record.description,
        length=len(sequence),
        gc_percent=percent,
        gc_bases=gc_bases,
        unambiguous_bases=unambiguous_bases,
        ambiguous_bases=ambiguous_bases,
    )
