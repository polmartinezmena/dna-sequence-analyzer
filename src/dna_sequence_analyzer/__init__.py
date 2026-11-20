"""DNA sequence analysis toolkit.

Stage 1 reads a FASTA file and reports, for every record, its id,
its length and its GC content. It can also estimate the melting
temperature of a short DNA oligo. Later stages (transcription, ORFs,
BLAST) are not part of this version.
"""

from dna_sequence_analyzer.gc import SequenceStats, analyze_fasta, gc_percent
from dna_sequence_analyzer.tm import (
    OLIGO_GUIDANCE_LENGTH,
    MeltingConditions,
    MeltingRecord,
    analyze_melting,
    melting_temperature,
)

__version__ = "0.1.1"

__all__ = [
    "MeltingConditions",
    "MeltingRecord",
    "OLIGO_GUIDANCE_LENGTH",
    "SequenceStats",
    "analyze_fasta",
    "analyze_melting",
    "gc_percent",
    "melting_temperature",
    "__version__",
]
