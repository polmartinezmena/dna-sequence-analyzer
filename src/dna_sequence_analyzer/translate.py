"""Transcription and six-frame translation using Biopython."""

from dataclasses import dataclass

from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq

FRAMES = ("+1", "+2", "+3", "-1", "-2", "-3")
TABLES = tuple(sorted(CodonTable.unambiguous_dna_by_id))


def _dna(sequence: str) -> Seq:
    dna = sequence.upper().replace("U", "T")
    invalid = set(dna).difference(IUPACData.ambiguous_dna_letters)
    if invalid:
        raise ValueError(f"Bases desconocidas: {', '.join(sorted(invalid))}")
    return Seq(dna)


def transcribe(sequence: str) -> str:
    return str(_dna(sequence).transcribe())


def reverse_complement(sequence: str) -> str:
    return str(_dna(sequence).reverse_complement())


@dataclass(frozen=True)
class FrameTranslation:
    frame: str
    protein: str

    @property
    def stops(self) -> int:
        return self.protein.count("*")

    @property
    def longest_stretch(self) -> int:
        return max(map(len, self.protein.split("*")))


def translate_frames(sequence: str, table: int = 1, frame: str = "all") -> list[FrameTranslation]:
    if table not in TABLES:
        raise ValueError(f"Tabla NCBI desconocida: {table}")
    if frame != "all" and frame not in FRAMES:
        raise ValueError(f"Marco desconocido: {frame}. Usa all o {', '.join(FRAMES)}.")
    forward = _dna(sequence)
    reverse = forward.reverse_complement()
    results = []
    for label in FRAMES if frame == "all" else (frame,):
        strand = forward if label.startswith("+") else reverse
        coding = strand[abs(int(label)) - 1:]
        complete = coding[:len(coding) // 3 * 3]
        results.append(FrameTranslation(label, str(complete.translate(table=table))))
    return results
