import warnings

import pytest
from typer.testing import CliRunner

from dna_sequence_analyzer.cli import app
from dna_sequence_analyzer.translate import reverse_complement, transcribe, translate_frames


def test_transcription_and_complement():
    assert transcribe("atgN") == "AUGN"
    assert reverse_complement("aRgt") == "ACYT"


def test_six_frames():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        frames = translate_frames("ATGAAATAG")
    assert [(f.frame, f.protein) for f in frames] == [
        ("+1", "MK*"), ("+2", "*N"), ("+3", "EI"),
        ("-1", "LFH"), ("-2", "YF"), ("-3", "IS"),
    ]


def test_tables_and_ambiguity():
    standard, = translate_frames("NNNTGAGCN", frame="+1")
    assert standard.protein == "X*A"
    mitochondrial, = translate_frames("NNNTGAGCN", table=2, frame="+1")
    assert mitochondrial.protein == "XWA"


@pytest.mark.parametrize("seq,stops,longest", [
    ("A", 0, 0), ("TAATGA", 2, 0), ("ATGTAAGCTGCTTAG", 2, 2),
])
def test_frame_stats(seq, stops, longest):
    item, = translate_frames(seq, frame="+1")
    assert (item.stops, item.longest_stretch) == (stops, longest)


@pytest.mark.parametrize("sequence,options", [("AT?", {}), ("ATG", {"table": 999}),
                                               ("ATG", {"frame": "0"})])
def test_invalid_input(sequence, options):
    with pytest.raises(ValueError):
        translate_frames(sequence, **options)


def test_cli_multiple_records(tmp_path):
    path = tmp_path / "records.fa"
    path.write_text(">uno\nTCACAT\n>dos\nCTA\n")
    result = CliRunner().invoke(app, ["translate", str(path), "--frame", "-1", "--table", "2"])
    assert result.exit_code == 0
    assert [line.split() for line in result.stdout.splitlines()[2:]] == [
        ["uno", "-1", "2", "0", "2"], ["dos", "-1", "0", "1", "0"],
    ]
