"""Slide a window along a sequence and summarize each piece.

A whole-record GC% hides local structure. A promoter, a CpG island or
the two halves of a bacterial chromosome can differ from the average.
These functions cut each FASTA record into windows of a fixed length,
moving forward by a fixed step, and score each window on its own.

Coordinates in the results are 1-based and inclusive, which is how a
sequence is usually numbered on paper: the first base is position 1.
A window that would run past the end of the record is not emitted.
If the record is shorter than the window, it contributes no rows.

CpG observed/expected uses the common ratio

    O/E = (number of CG dinucleotides * window length) / (number of C * number of G)

Only a literal C followed by a literal G counts as CpG. S is not
treated as both, because that would invent dinucleotides. When the
window has no C or no G, the expected count is zero and O/E is
undefined.

GC skew is (G - C) / (G + C), again with literal G and C only. It runs
from -1 (C and no G) to +1 (G and no C). The cumulative skew is the
running sum of those window values, restarted for every record. In
many bacterial chromosomes the leading strand is richer in G than in
C, so the cumulative curve falls along one replichore and rises along
the other. The minimum is often near the replication origin and the
maximum near the terminus. That pattern is a clue, not a proof, and
it is a bacterial-chromosome tool rather than a primer tool.
"""

from dataclasses import dataclass
from pathlib import Path

from dna_sequence_analyzer.gc import _count_bases, gc_percent, read_fasta


@dataclass(frozen=True)
class GCWindow:
    """GC% and CpG observed/expected for one window of one record.

    Attributes:
        id: FASTA identifier of the parent record.
        start: First base of the window, 1-based.
        end: Last base of the window, 1-based and inclusive.
        gc_percent: GC% of this window, or ``None`` when it has no
            countable bases.
        cpg_oe: CpG observed/expected, or ``None`` when the window has
            no C or no G.
        ambiguous_bases: IUPAC ambiguity codes inside the window.
    """

    id: str
    start: int
    end: int
    gc_percent: float | None
    cpg_oe: float | None
    ambiguous_bases: int


@dataclass(frozen=True)
class SkewWindow:
    """GC skew of one window, plus the cumulative sum for its record.

    Attributes:
        id: FASTA identifier of the parent record.
        start: First base of the window, 1-based.
        end: Last base of the window, 1-based and inclusive.
        gc_skew: ``(G - C) / (G + C)``, or ``None`` when the window
            has neither G nor C.
        cumulative_skew: Sum of the defined window skews from the
            start of this record through this window. ``None`` until
            the first window that has a defined skew.
    """

    id: str
    start: int
    end: int
    gc_skew: float | None
    cumulative_skew: float | None


def gc_windows(path: str | Path, window: int, step: int) -> list[GCWindow]:
    """Score every complete window of every FASTA record.

    Args:
        path: FASTA file.
        window: Window length in bases. At least 1.
        step: Distance from one window start to the next. At least 1.
            Equal to ``window`` means the windows do not overlap.
            Smaller than ``window`` means they overlap.

    Raises:
        FileNotFoundError: ``path`` does not exist.
        ValueError: ``path`` is not a FASTA file, or ``window`` or
            ``step`` is below 1.
    """
    return gc_windows_for_records(read_fasta(path), window, step)


def gc_windows_for_records(records, window: int, step: int) -> list[GCWindow]:
    """Same score as :func:`gc_windows`, for records already in memory.

    The web app uses this so it does not write a temporary FASTA file
    and does not reimplement the window math.
    """
    _check_span(window, step)
    rows: list[GCWindow] = []
    for record in records:
        sequence = str(record.seq)
        for start, end in _window_bounds(len(sequence), window, step):
            chunk = sequence[start:end]
            _gc_bases, _at_bases, ambiguous = _count_bases(chunk)
            rows.append(
                GCWindow(
                    id=record.id,
                    start=start + 1,
                    end=end,
                    gc_percent=gc_percent(chunk),
                    cpg_oe=cpg_observed_expected(chunk),
                    ambiguous_bases=ambiguous,
                )
            )
    return rows


def skew_windows(path: str | Path, window: int, step: int) -> list[SkewWindow]:
    """Compute window GC skew and cumulative skew for every record.

    The cumulative sum restarts at each FASTA record. Windows with no
    G and no C do not move the cumulative sum.

    Raises:
        FileNotFoundError: ``path`` does not exist.
        ValueError: ``path`` is not a FASTA file, or ``window`` or
            ``step`` is below 1.
    """
    return skew_windows_for_records(read_fasta(path), window, step)


def skew_windows_for_records(records, window: int, step: int) -> list[SkewWindow]:
    """Same score as :func:`skew_windows`, for records already in memory."""
    _check_span(window, step)
    rows: list[SkewWindow] = []
    for record in records:
        sequence = str(record.seq)
        cumulative = 0.0
        seen_defined = False
        for start, end in _window_bounds(len(sequence), window, step):
            skew = gc_skew(sequence[start:end])
            if skew is not None:
                cumulative += skew
                seen_defined = True
            rows.append(
                SkewWindow(
                    id=record.id,
                    start=start + 1,
                    end=end,
                    gc_skew=skew,
                    cumulative_skew=cumulative if seen_defined else None,
                )
            )
    return rows


def skew_extrema(rows: list[SkewWindow]) -> tuple[SkewWindow | None, SkewWindow | None]:
    """Return the windows with the lowest and highest cumulative skew.

    Both are ``None`` when no window has a defined cumulative value.
    If several windows share the minimum or the maximum, the first one
    in the list is returned. On a bacterial chromosome those two points
    are the usual candidates for the replication origin (minimum) and
    the terminus (maximum). They are candidates, not a proof.
    """
    defined = [row for row in rows if row.cumulative_skew is not None]
    if not defined:
        return None, None
    lowest = min(defined, key=lambda row: row.cumulative_skew)
    highest = max(defined, key=lambda row: row.cumulative_skew)
    return lowest, highest


def cpg_observed_expected(sequence: str) -> float | None:
    """Return CpG observed/expected, or ``None`` if it is undefined.

    ``O/E = (CG dinucleotides * length) / (C count * G count)``.

    Counts are case-insensitive. Only the dinucleotide CG counts, and
    only the letters C and G enter the expected value. S is left out.
    """
    seq = "".join(sequence.split()).upper()
    observed = sum(seq[index : index + 2] == "CG" for index in range(len(seq) - 1))
    cytosines = seq.count("C")
    guanines = seq.count("G")
    if cytosines == 0 or guanines == 0 or not seq:
        return None
    return (observed * len(seq)) / (cytosines * guanines)


def gc_skew(sequence: str) -> float | None:
    """Return ``(G - C) / (G + C)``, or ``None`` when both counts are zero.

    Literal G and C only. S is excluded: it is G or C, but the sign of
    the skew depends on which one, and that is unknown.
    """
    seq = "".join(sequence.split()).upper()
    guanines = seq.count("G")
    cytosines = seq.count("C")
    if guanines + cytosines == 0:
        return None
    return (guanines - cytosines) / (guanines + cytosines)


def _check_span(window: int, step: int) -> None:
    for name, value in (("window", window), ("step", step)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be an integer of at least 1.")


def _window_bounds(length: int, window: int, step: int):
    """Yield ``(start, end)`` half-open slices that fit inside ``length``."""
    start = 0
    while start + window <= length:
        yield start, start + window
        start += step
