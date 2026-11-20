"""Command-line interface for the DNA sequence analyzer.

The only command in this version is ``gc``, which prints a table of
sequence id, length and GC content. Later stages will add more commands.
"""

import csv
import sys
from pathlib import Path
from typing import TextIO

import typer

from dna_sequence_analyzer import __version__
from dna_sequence_analyzer.gc import SequenceStats, analyze_fasta

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "DNA sequence analysis toolkit. "
        "Stage 1 reads a FASTA file and reports length and GC content."
    ),
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"dna-analyzer {__version__}")
        raise typer.Exit()


@app.callback()
def callback(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the program version and exit.",
    ),
) -> None:
    """Analyze DNA sequences from the command line.

    ``--version`` is handled before this function body runs. The
    parameter is what makes Typer show the flag in the help text.
    """


@app.command("gc")
def gc(
    fasta_file: Path = typer.Argument(
        ...,
        metavar="FASTA",
        help="FASTA file with one or more DNA (or RNA) records.",
    ),
    as_csv: bool = typer.Option(
        False,
        "--csv",
        help="Print id,length,gc_percent as CSV instead of a text table.",
    ),
) -> None:
    """Print the length and GC content of each FASTA record.

    GC% counts G, C and S as GC, and A, T, U and W as known bases that
    are not GC. Other ambiguity codes, such as N, are excluded. Length
    still includes every symbol in the sequence. GC% is shown as NA
    when a record has no bases that can be counted.
    """
    try:
        stats = analyze_fasta(fasta_file)
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    if as_csv:
        _write_csv(stats, sys.stdout)
    else:
        typer.echo(_format_table(stats))

    if any(item.gc_percent is None for item in stats):
        typer.echo(
            "Note: GC% is NA when a record has no bases that can be counted "
            "(A, C, G, T, U, S or W).",
            err=True,
        )


def main() -> None:
    """Entry point used by the ``dna-analyzer`` console script."""
    app()


def _format_gc(gc_percent: float | None) -> str:
    if gc_percent is None:
        return "NA"
    return f"{gc_percent:.2f}"


def _format_table(stats: list[SequenceStats]) -> str:
    """Render records as an aligned text table.

    Identifiers are left-aligned. Length and GC% are right-aligned so
    the digits line up when the table is read down a column.
    """
    header = ("id", "length", "gc_percent")
    rows = [
        (item.id, str(item.length), _format_gc(item.gc_percent)) for item in stats
    ]
    widths = [
        max(len(header[column]), *(len(row[column]) for row in rows))
        for column in range(3)
    ]

    def render(columns: tuple[str, str, str]) -> str:
        id_text, length_text, gc_text = columns
        return (
            f"{id_text.ljust(widths[0])}  "
            f"{length_text.rjust(widths[1])}  "
            f"{gc_text.rjust(widths[2])}"
        )

    separator = "  ".join("-" * width for width in widths)
    lines = [render(header), separator]
    lines.extend(render(row) for row in rows)
    return "\n".join(lines)


def _write_csv(stats: list[SequenceStats], destination: TextIO) -> None:
    """Write records as CSV with a header row.

    ``gc_percent`` is empty when it is undefined, so a spreadsheet does
    not treat the text NA as a number. The text table uses NA instead,
    where an empty cell would be easy to miss.
    """
    writer = csv.writer(destination, lineterminator="\n")
    writer.writerow(["id", "length", "gc_percent"])
    for item in stats:
        gc_text = "" if item.gc_percent is None else f"{item.gc_percent:.2f}"
        writer.writerow([item.id, item.length, gc_text])
