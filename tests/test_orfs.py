import pytest

from dna_sequence_analyzer.orfs import find_orfs
from dna_sequence_analyzer.cli import app
from typer.testing import CliRunner


@pytest.mark.parametrize(
    ("sequence", "frame", "start", "end"),
    [("catgaaataagg", "+2", 2, 10), ("ccttatttcatg", "-2", 3, 11)],
)
def test_coordinates_reference_original_sequence(sequence, frame, start, end):
    hits = find_orfs(sequence, min_aa=2)
    assert [(x.frame, x.start, x.end, x.protein) for x in hits] == [
        (frame, start, end, "MK")
    ]


def test_alternative_start_and_nested_hits():
    sequence = "GTGATGNNNTAA"
    assert [x.protein for x in find_orfs(sequence, table=11, min_aa=1, alternative_starts=False)] == ["MX"]
    hits = find_orfs(sequence, table=11, min_aa=1, alternative_starts=True, include_nested=True)
    assert [(x.start_codon, x.protein) for x in hits] == [("GTG", "MMX"), ("ATG", "MX")]
    hits = find_orfs(sequence, table=11, min_aa=1, alternative_starts=True, include_nested=False)
    assert [x.protein for x in hits] == ["MMX"]


@pytest.mark.parametrize("sequence", ["", "NNN", "ATGAAA", "ATGAAATAA"])
def test_missing_stop_or_below_threshold(sequence):
    assert find_orfs(sequence, min_aa=3) == []


def test_rejects_invalid_bases():
    with pytest.raises(ValueError):
        find_orfs("ATG?")


def test_open_end_on_both_strands():
    from Bio.Seq import Seq

    for sequence, frame, start, end in [("ATGAAAC", "+1", 1, 6), ("GTTTCAT", "-1", 2, 7)]:
        hits = find_orfs(sequence, allow_no_stop=True, min_aa=2)
        assert len(hits) == 1
        hit = hits[0]
        coding = Seq(sequence[start - 1:end])
        if frame == "-1":
            coding = coding.reverse_complement()
        assert (hit.frame, hit.start, hit.end, hit.has_stop) == (frame, start, end, False)
        assert hit.protein == str(coding.translate())


@pytest.mark.parametrize(("kind", "expected"), [
    ("csv", "demo,-2,3,11,2,MK"),
    ("fasta", ">orf1 demo -2 3-11\nMK"),
    ("gff", "demo\tdna-analyzer\topen_reading_frame\t3\t11\t.\t-\t.\tID=orf1"),
])
def test_cli_exports_reverse_orf(tmp_path, kind, expected):
    fasta = tmp_path / "demo.fasta"
    fasta.write_text(">demo\nCCTTATTTCATG\n", encoding="ascii")
    result = CliRunner().invoke(app, ["orfs", str(fasta), "--min-aa", "2", "--format", kind])
    assert result.exit_code == 0, result.output
    assert expected in result.stdout


def test_insulin_default_starts():
    from dna_sequence_analyzer.gc import read_fasta
    from dna_sequence_analyzer.web_input import bundled_example_fasta

    insulin = read_fasta(bundled_example_fasta())[0]
    hits = find_orfs(str(insulin.seq), min_aa=30)
    assert [(h.frame, h.start, h.end, len(h.protein)) for h in hits] == [
        ("+3", 60, 392, 110), ("-3", 206, 355, 49),
    ]


def test_seeded_orfs_match_coordinates_on_both_strands():
    from random import Random

    from Bio.Seq import Seq

    sequence = Seq("".join(Random(42).choices("ACGT", k=6000)))
    hits = find_orfs(str(sequence))
    assert {h.strand for h in hits} == {"+", "-"}
    for hit in hits:
        fragment = sequence[hit.start - 1:hit.end]
        if hit.strand == "-":
            fragment = fragment.reverse_complement()
        assert len(fragment) % 3 == 0
        assert fragment[:3] == hit.start_codon == "ATG"
        assert str(fragment.translate()) == hit.protein + "*"
        assert hit.start_codon_position == (hit.start if hit.strand == "+" else hit.end)


def test_cli_alternative_start_and_codon_columns(tmp_path):
    import csv

    fasta = tmp_path / "alternativo.fasta"
    fasta.write_text(">demo\nGTGAAATAA\n", encoding="ascii")
    args = ["orfs", str(fasta), "--table", "11", "--min-aa", "2"]
    runner = CliRunner()
    assert not list(csv.DictReader(runner.invoke(app, args).stdout.splitlines()))
    result = runner.invoke(app, [*args, "--alt-starts"])
    assert result.exit_code == 0, result.output
    row, = csv.DictReader(result.stdout.splitlines())
    assert [row[key] for key in ("Hebra", "Codón de inicio", "Posición del codón de inicio")] == ["+", "GTG", "1"]
    assert row["protein"] == "MK"
