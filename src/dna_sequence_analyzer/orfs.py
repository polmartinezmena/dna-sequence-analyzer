"""Find complete start-to-stop ORFs in linear DNA."""

import csv
from dataclasses import dataclass
from io import StringIO
from urllib.parse import quote

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord

from .translate import FRAMES, _dna


@dataclass(frozen=True)
class ORF:
    """Coordinates are 1-based and inclusive, including the stop when present."""

    frame: str
    start: int
    end: int
    protein: str
    has_stop: bool = True
    start_codon: str = "ATG"

    @property
    def strand(self) -> str:
        return self.frame[0]

    @property
    def start_codon_position(self) -> int:
        return self.start if self.strand == "+" else self.end


def find_orfs(
    sequence: str,
    *,
    table: int = 1,
    min_aa: int = 1,
    alternative_starts: bool = False,
    include_nested: bool = False,
    allow_no_stop: bool = False,
) -> list[ORF]:
    """Find candidates in all six frames using Biopython genetic codes."""
    if min_aa < 1:
        raise ValueError("La longitud mínima debe ser al menos 1 aminoácido.")
    try:
        code = CodonTable.unambiguous_dna_by_id[table]
    except KeyError:
        raise ValueError(f"Tabla NCBI desconocida: {table}") from None
    stops = set(code.stop_codons).difference(code.forward_table)
    if not stops:
        raise ValueError("La tabla no tiene codones de parada inequívocos.")
    starts = set(code.start_codons) if alternative_starts else {"ATG"}
    forward = _dna(sequence)
    reverse = forward.reverse_complement()
    size = len(forward)
    hits = []
    for frame in FRAMES:
        positive = frame.startswith("+")
        strand = forward if positive else reverse
        pending = []
        offset = abs(int(frame)) - 1
        end_of_frame = offset + max(0, (size - offset) // 3) * 3
        for pos in range(offset, end_of_frame + int(allow_no_stop), 3):
            has_stop = pos < end_of_frame
            codon = str(strand[pos:pos + 3]) if has_stop else None
            if not has_stop or codon in stops:
                for start in pending:
                    if pos - start < 3 * min_aa:
                        break
                    end = pos + 3 if has_stop else pos
                    protein = "M" + str(strand[start + 3:pos].translate(table=table))
                    left, right = (start, end) if positive else (size - end, size - start)
                    hits.append(ORF(frame, left + 1, right, protein, has_stop, str(strand[start:start + 3])))
                pending.clear()
            elif codon in starts and (include_nested or not pending):
                pending.append(pos)
    return hits


def export_orfs(records: list[tuple[str, ORF]], kind: str) -> str:
    """Serialize results using one-based inclusive coordinates."""
    output = StringIO()
    if kind == "fasta":
        write(
            (SeqRecord(Seq(hit.protein), id=f"orf{i}", description=f"{name} {hit.frame} {hit.start}-{hit.end}")
             for i, (name, hit) in enumerate(records, 1)),
            output, "fasta",
        )
    elif kind == "gff":
        output.write("##gff-version 3\n")
        for i, (name, hit) in enumerate(records, 1):
            seqid = quote(name, safe=".-_:")
            output.write(
                f"{seqid}\tdna-analyzer\topen_reading_frame\t{hit.start}\t{hit.end}"
                f"\t.\t{hit.frame[0]}\t.\tID=orf{i}\n"
            )
    elif kind == "csv":
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("id", "frame", "start_1based", "end_inclusive", "length_aa", "protein", "has_stop",
                         "Hebra", "Cod?n de inicio", "Posici?n del cod?n de inicio"))
        writer.writerows((name, hit.frame, hit.start, hit.end, len(hit.protein), hit.protein, hit.has_stop, hit.strand.replace("-", "?"),
                          hit.start_codon, hit.start_codon_position)
                         for name, hit in records)
    else:
        raise ValueError(f"Formato desconocido: {kind}")
    return output.getvalue()
