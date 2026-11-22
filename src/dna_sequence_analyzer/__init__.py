"""DNA sequence analysis toolkit.

The package reads FASTA, reports length and GC content, slides a
window for local GC and GC skew, and estimates oligo melting
temperature. The toolkit also transcribes DNA and translates all six frames.
"""

from dna_sequence_analyzer.gc import (
    GC_HIGH_PERCENT,
    GC_LOW_PERCENT,
    SequenceStats,
    analyze_fasta,
    gc_percent,
    nucleic_acid_kind,
)
from dna_sequence_analyzer.tm import (
    OLIGO_GUIDANCE_LENGTH,
    MeltingConditions,
    MeltingRecord,
    analyze_melting,
    melting_temperature,
)
from dna_sequence_analyzer.windows import (
    GCWindow,
    SkewWindow,
    gc_windows,
    skew_windows,
)

__version__ = "1.0.0"

__all__ = [
    "GCWindow",
    "GC_HIGH_PERCENT",
    "GC_LOW_PERCENT",
    "MeltingConditions",
    "MeltingRecord",
    "OLIGO_GUIDANCE_LENGTH",
    "SequenceStats",
    "SkewWindow",
    "analyze_fasta",
    "analyze_melting",
    "gc_percent",
    "gc_windows",
    "melting_temperature",
    "nucleic_acid_kind",
    "skew_windows",
    "__version__",
]
