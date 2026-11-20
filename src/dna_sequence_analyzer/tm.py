"""Melting temperature of a short DNA oligo.

The melting temperature (Tm) is the temperature at which half of the
duplexes in a solution have separated into single strands. It is not a
property of the GC percentage. Two oligos can have the same GC% and
different Tm values, because stability depends on the order of the
bases.

The physical reason is base stacking. Along a duplex, each base pair
sits on top of the next one, and that contact is the main contribution
to duplex stability. A GC-containing step is often a stronger stack
than an AT-only step, but "often" is not "always", and the particular
neighbors matter. Hydrogen bonds are real (three between G and C, two
between A and T), yet they contribute less to the *net* stability than
the bond count suggests: a base that is not paired with its partner is
still hydrogen-bonded to water. Opening the helix trades base-base
bonds for base-water bonds.

Nearest-neighbor models take that order dependence literally. They add
a measured enthalpy and entropy for every dinucleotide step (AA, AT,
GC, CG, and so on), then adjust for salt and for how much oligo is in
the tube. This module calls Biopython's ``Bio.SeqUtils.MeltingTemp.Tm_NN``
with the DNA/DNA table DNA_NN3 (Allawi and SantaLucia, 1997), which is
Biopython's default.

The calculation is meant for short DNA oligos and PCR primers, roughly
60 nucleotides or fewer. It is not a melting temperature for a long
mRNA or a genome. A long sequence still produces a number, and the
command-line tool says so, but that number should not be read as a
PCR condition.
"""

from dataclasses import dataclass
from pathlib import Path

from Bio.SeqUtils.MeltingTemp import Tm_NN

from dna_sequence_analyzer.gc import read_fasta

# Past this length the nearest-neighbor oligo model is the wrong tool.
# Typical PCR primers are about 18-25 nt. Sixty is a generous ceiling
# for "short oligo", not a biological phase transition.
OLIGO_GUIDANCE_LENGTH = 60

# DNA bases, plus U so an RNA-style letter can be read as T.
# Anything else, including N, is rejected. Biopython's own checker
# would silently delete those letters before calculating, which would
# report the Tm of a different molecule.
_ACCEPTED_BASES = frozenset("ACGTU")


@dataclass(frozen=True)
class MeltingConditions:
    """Solution conditions passed through to ``Tm_NN``.

    The defaults are Biopython's defaults. Ion concentrations are
    millimolar. Strand concentrations are nanomolar.

    Attributes:
        na_mM: Na+ concentration. Default 50 mM. Cations gather around
            the negatively charged phosphate backbone and screen the
            repulsion between the two strands, which raises Tm. Mg2+
            does this more strongly, per ion, than Na+.
        k_mM: K+ concentration. Default 0.
        tris_mM: Tris buffer concentration. Default 0.
        mg_mM: Mg2+ concentration. Default 0.
        dntps_mM: dNTP concentration. Default 0. In a PCR mix, free
            magnesium is what is left after dNTPs bind some of it.
            This calculator does not guess a PCR recipe; pass the
            concentrations you want ``Tm_NN`` to use.
        dnac1_nM: Concentration of the more abundant strand, typically
            the primer. Default 25 nM.
        dnac2_nM: Concentration of the less abundant strand, typically
            the template. Default 25 nM. Biopython's note on Primer3Plus:
            a total oligo concentration of 50 nM is ``dnac1=25`` and
            ``dnac2=25``.
        saltcorr: Which salt correction ``Tm_NN`` applies. Default 5,
            the SantaLucia 1998 correction to the entropy. Biopython
            converts the millimolar ion concentrations to molar before
            that formula. ``0`` turns the correction off. Allowed values
            are 0 through 7; they are the methods documented on
            ``Bio.SeqUtils.MeltingTemp.salt_correction``.
        selfcomp: When True, the oligo is treated as binding to itself
            and ``dnac2_nM`` is ignored. Default False.
    """

    na_mM: float = 50
    k_mM: float = 0
    tris_mM: float = 0
    mg_mM: float = 0
    dntps_mM: float = 0
    dnac1_nM: float = 25
    dnac2_nM: float = 25
    saltcorr: int = 5
    selfcomp: bool = False

    def __post_init__(self) -> None:
        for name in ("na_mM", "k_mM", "tris_mM", "mg_mM", "dntps_mM", "dnac1_nM", "dnac2_nM"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} cannot be negative.")
        if isinstance(self.saltcorr, bool) or not isinstance(self.saltcorr, int):
            raise ValueError("saltcorr must be an integer from 0 to 7.")
        if self.saltcorr not in range(8):
            raise ValueError("saltcorr must be an integer from 0 to 7.")
        if self.selfcomp:
            if self.dnac1_nM <= 0:
                raise ValueError("dnac1_nM must be positive.")
        elif self.dnac1_nM - (self.dnac2_nM / 2.0) <= 0:
            raise ValueError(
                "dnac1_nM must be greater than dnac2_nM / 2. "
                "Tm_NN takes the logarithm of that concentration. "
                "For a self-complementary oligo, set selfcomp=True."
            )


@dataclass(frozen=True)
class MeltingRecord:
    """Nearest-neighbor Tm of one sequence.

    Attributes:
        id: Sequence name. For a FASTA record this is the first word
            of the header. For a sequence typed on the command line it
            is ``sequence``.
        description: Full FASTA header, without the leading ``>``.
            Empty when the sequence did not come from a file.
        length: Number of symbols in the sequence as supplied, after
            removing whitespace. Ambiguous letters are included.
        tm_celsius: Melting temperature in degrees Celsius, or ``None``
            when the sequence cannot be evaluated honestly.
        sequence_is_long: True when ``length`` is above
            :data:`OLIGO_GUIDANCE_LENGTH`. The Tm is still reported.
            A long value is not a primer melting temperature.
        contains_uracil: True when the sequence contains U. U is read
            as T for the DNA/DNA table. The number is not an RNA
            folding temperature, and G-U wobble pairs are not modeled.
    """

    id: str
    description: str
    length: int
    tm_celsius: float | None
    sequence_is_long: bool
    contains_uracil: bool


def melting_temperature(
    sequence: str,
    conditions: MeltingConditions | None = None,
) -> float:
    """Return the nearest-neighbor Tm of one DNA oligo, in Celsius.

    Whitespace is ignored and case does not matter. ``U`` is read as
    ``T``, so a sequence written with RNA letters can be scored with
    the DNA/DNA table. That is a convenience for mixed notation, not
    an RNA melting model: the RNA thermodynamic tables are not used.

    Args:
        sequence: Oligo sequence using A, C, G, T or U only.
        conditions: Salt, strand concentration and salt-correction
            method. ``None`` uses :class:`MeltingConditions` defaults,
            which match ``Tm_NN``.

    Returns:
        Melting temperature in degrees Celsius.

    Raises:
        ValueError: The sequence is shorter than 2 bases, contains a
            letter other than A, C, G, T or U, or ``conditions`` is not
            usable. Ambiguous bases are not stripped. Dropping an ``N``
            would silently calculate a different oligo.
    """
    settings = conditions if conditions is not None else MeltingConditions()
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise ValueError("Sequence is empty. Nearest-neighbor Tm needs at least 2 bases.")
    unexpected = sorted({base for base in cleaned if base not in _ACCEPTED_BASES})
    if unexpected:
        letters = ", ".join(unexpected)
        raise ValueError(
            "Tm_NN needs an unambiguous DNA sequence (A, C, G, T or U). "
            f"Found {letters}. Ambiguous letters are not removed, because "
            "dropping them would score a different oligo."
        )
    if len(cleaned) < 2:
        raise ValueError("Nearest-neighbor Tm needs at least 2 bases.")
    # DNA_NN3 is a DNA/DNA table. U is uracil; the DNA equivalent is T.
    dna = cleaned.replace("U", "T")
    return float(
        Tm_NN(
            dna,
            dnac1=settings.dnac1_nM,
            dnac2=settings.dnac2_nM,
            selfcomp=settings.selfcomp,
            Na=settings.na_mM,
            K=settings.k_mM,
            Tris=settings.tris_mM,
            Mg=settings.mg_mM,
            dNTPs=settings.dntps_mM,
            saltcorr=settings.saltcorr,
        )
    )


def analyze_melting(
    path: str | Path,
    conditions: MeltingConditions | None = None,
) -> list[MeltingRecord]:
    """Calculate a nearest-neighbor Tm for every record in a FASTA file.

    A record that is empty, shorter than 2 bases, or contains letters
    other than A, C, G, T or U is kept in the result with
    ``tm_celsius is None``. The other records are still calculated.

    Args:
        path: Path to a FASTA file.
        conditions: Passed to :func:`melting_temperature`. ``None`` uses
            the Biopython defaults.

    Raises:
        FileNotFoundError: ``path`` does not exist.
        ValueError: ``path`` is not a readable FASTA file, or
            ``conditions`` is not usable.
    """
    settings = conditions if conditions is not None else MeltingConditions()
    records = []
    for record in read_fasta(path):
        sequence = str(record.seq)
        records.append(_record_from_sequence(record.id, record.description, sequence, settings))
    return records


def melting_record_for_sequence(
    sequence: str,
    conditions: MeltingConditions | None = None,
    *,
    sequence_id: str = "sequence",
) -> MeltingRecord:
    """Calculate the Tm of one sequence that did not come from a file."""
    settings = conditions if conditions is not None else MeltingConditions()
    return _record_from_sequence(sequence_id, "", sequence, settings)


def _record_from_sequence(
    sequence_id: str,
    description: str,
    sequence: str,
    conditions: MeltingConditions,
) -> MeltingRecord:
    cleaned = "".join(sequence.split())
    try:
        tm = melting_temperature(cleaned, conditions)
    except ValueError:
        tm = None
    return MeltingRecord(
        id=sequence_id,
        description=description,
        length=len(cleaned),
        tm_celsius=tm,
        sequence_is_long=len(cleaned) > OLIGO_GUIDANCE_LENGTH,
        contains_uracil="U" in cleaned.upper(),
    )
