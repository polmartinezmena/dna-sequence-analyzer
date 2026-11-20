"""DNA sequence analysis toolkit.

Stage 1 reads a FASTA file and reports, for every record, its id,
its length and its GC content. Later stages (transcription, ORFs,
BLAST) are not part of this version.
"""

from dna_sequence_analyzer.gc import SequenceStats, analyze_fasta, gc_percent

__version__ = "0.1.0"

__all__ = [
    "SequenceStats",
    "analyze_fasta",
    "gc_percent",
    "__version__",
]
