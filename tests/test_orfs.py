import pytest

from dna_sequence_analyzer.orfs import find_orfs
from dna_sequence_analyzer.cli import app
from typer.testing import CliRunner


@pytest.mark.parametrize(
    ("sequence", "frame", "start", "end"),
    [("catgaaataagg", "+2", 1, 10), ("ccttatttcatg", "-2", 2, 11)],
)
def test_coordinates_reference_original_sequence(sequence, frame, start, end):
    hits = find_orfs(sequence, min_aa=2)
    assert [(x.frame, x.start, x.end, x.protein) for x in hits] == [
        (frame, start, end, "MK")
    ]


def test_alternative_start_and_nested_hits():
    sequence = "GTGATGNNNTAA"
    assert [x.protein for x in find_orfs(sequence, table=11, min_aa=1)] == ["MX"]
    hits = find_orfs(sequence, table=11, min_aa=1, alternative_starts=True)
    assert [x.protein for x in hits] == ["MMX", "MX"]
    hits = find_orfs(sequence, table=11, min_aa=1, alternative_starts=True, include_nested=False)
    assert [x.protein for x in hits] == ["MMX"]


@pytest.mark.parametrize("sequence", ["", "NNN", "ATGAAA", "ATGAAATAA"])
def test_missing_stop_or_below_threshold(sequence):
    assert find_orfs(sequence, min_aa=3) == []


def test_rejects_invalid_bases():
    with pytest.raises(ValueError):
        find_orfs("ATG?")


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
