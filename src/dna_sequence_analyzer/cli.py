"""Command-line interface for the DNA sequence analyzer.

``gc`` prints sequence id, length and GC content. ``tm`` prints a
nearest-neighbor melting temperature for short DNA oligos. Later stages
will add more commands.
"""

import csv
import sys
from pathlib import Path
from typing import TextIO

import typer

from dna_sequence_analyzer import __version__
from dna_sequence_analyzer.gc import analyze_fasta
from dna_sequence_analyzer.tm import (
    OLIGO_GUIDANCE_LENGTH,
    MeltingConditions,
    MeltingRecord,
    analyze_melting,
    melting_record_for_sequence,
)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "DNA sequence analysis toolkit. "
        "Reads FASTA and reports length, GC content, and oligo melting temperature."
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

    header = ("id", "length", "gc_percent")
    rows = [
        (
            item.id,
            str(item.length),
            "" if item.gc_percent is None else f"{item.gc_percent:.2f}",
        )
        for item in stats
    ]
    if as_csv:
        _write_csv(header, rows, sys.stdout)
    else:
        display_rows = [
            (item_id, length, gc_text or "NA") for item_id, length, gc_text in rows
        ]
        typer.echo(_format_table(header, display_rows))

    if any(item.gc_percent is None for item in stats):
        typer.echo(
            "Note: GC% is NA when a record has no bases that can be counted "
            "(A, C, G, T, U, S or W).",
            err=True,
        )


@app.command("tm")
def tm(
    fasta_file: Path | None = typer.Argument(
        None,
        metavar="[FASTA]",
        help="FASTA file of oligos. Omit it when you pass --sequence.",
    ),
    sequence: str | None = typer.Option(
        None,
        "--sequence",
        "-s",
        help="One oligo (A, C, G, T or U) instead of a FASTA file.",
    ),
    na: float = typer.Option(50, "--na", help="Na+ concentration in mM."),
    k: float = typer.Option(0, "--k", help="K+ concentration in mM."),
    tris: float = typer.Option(0, "--tris", help="Tris concentration in mM."),
    mg: float = typer.Option(0, "--mg", help="Mg2+ concentration in mM."),
    dntps: float = typer.Option(0, "--dntps", help="dNTP concentration in mM."),
    dnac1: float = typer.Option(
        25,
        "--dnac1",
        help="Concentration of the more abundant strand, usually the primer, in nM.",
    ),
    dnac2: float = typer.Option(
        25,
        "--dnac2",
        help="Concentration of the less abundant strand, usually the template, in nM.",
    ),
    saltcorr: int = typer.Option(
        5,
        "--saltcorr",
        help="Salt correction passed to Tm_NN (0-7). 5 is SantaLucia 1998. 0 disables it.",
        min=0,
        max=7,
    ),
    self_complementary: bool = typer.Option(
        False,
        "--self-complementary",
        help="Treat the oligo as binding to itself. dnac2 is then ignored.",
    ),
    as_csv: bool = typer.Option(
        False,
        "--csv",
        help="Print id,length,tm_celsius as CSV instead of a text table.",
    ),
) -> None:
    """Estimate the melting temperature of short DNA oligos.

    Uses Biopython's nearest-neighbor function Tm_NN and the DNA/DNA
    table DNA_NN3 (Allawi and SantaLucia, 1997). Defaults match Tm_NN:
    50 mM Na+, no K+, Tris, Mg2+ or dNTPs, 25 nM of each strand, and
    salt correction 5.

    This is for oligos and primers of about 60 nt or fewer. A long
    mRNA still gets a number, but that number is not a PCR melting
    temperature. The Tm depends on sequence order, not only on GC%.
    """
    if fasta_file is None and sequence is None:
        typer.echo("Error: pass a FASTA file or --sequence.", err=True)
        raise typer.Exit(code=1)
    if fasta_file is not None and sequence is not None:
        typer.echo("Error: pass a FASTA file or --sequence, not both.", err=True)
        raise typer.Exit(code=1)
    try:
        conditions = MeltingConditions(
            na_mM=na,
            k_mM=k,
            tris_mM=tris,
            mg_mM=mg,
            dntps_mM=dntps,
            dnac1_nM=dnac1,
            dnac2_nM=dnac2,
            saltcorr=saltcorr,
            selfcomp=self_complementary,
        )
        if sequence is not None:
            records = [melting_record_for_sequence(sequence, conditions)]
        else:
            assert fasta_file is not None
            records = analyze_melting(fasta_file, conditions)
    except FileNotFoundError as exc:
        typer.echo(
            f"Error: {exc}. If this was a sequence rather than a path, pass it with --sequence.",
            err=True,
        )
        raise typer.Exit(code=1) from exc
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    header = ("id", "length", "tm_celsius")
    rows = [
        (item.id, str(item.length), _format_optional(item.tm_celsius)) for item in records
    ]
    if as_csv:
        _write_csv(header, [_csv_tm_row(item) for item in records], sys.stdout)
    else:
        typer.echo(_format_table(header, rows))
    if any(item.tm_celsius is None for item in records):
        typer.echo(
            "Note: Tm is NA when a sequence is shorter than 2 bases or contains "
            "letters other than A, C, G, T or U. Ambiguous bases are not removed.",
            err=True,
        )
    if any(item.sequence_is_long for item in records):
        typer.echo(
            "Note: Tm_NN is meant for short DNA oligos and primers "
            f"(about {OLIGO_GUIDANCE_LENGTH} nt or fewer), not for long molecules "
            "such as mRNA. A longer record is still calculated, but do not use "
            "that number as a PCR melting temperature.",
            err=True,
        )


def main() -> None:
    """Entry point used by the ``dna-analyzer`` console script."""
    app()


def _format_optional(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.2f}"


def _format_table(header: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    """Render an aligned text table.

    The first column is left-aligned. The other columns are
    right-aligned so the digits line up when the table is read down
    a column.
    """
    widths = [
        max(len(header[column]), *(len(row[column]) for row in rows))
        for column in range(len(header))
    ]

    def render(columns: tuple[str, ...]) -> str:
        cells = [
            text.ljust(widths[index]) if index == 0 else text.rjust(widths[index])
            for index, text in enumerate(columns)
        ]
        return "  ".join(cells)

    separator = "  ".join("-" * width for width in widths)
    lines = [render(header), separator]
    lines.extend(render(row) for row in rows)
    return "\n".join(lines)


def _write_csv(
    header: tuple[str, ...],
    rows: list[tuple[str, ...]],
    destination: TextIO,
) -> None:
    """Write a header and rows as CSV.

    Undefined numbers are empty cells, so a spreadsheet does not treat
    the text NA as a number. The text table uses NA instead, where an
    empty cell would be easy to miss.
    """
    writer = csv.writer(destination, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)


def _csv_tm_row(item: MeltingRecord) -> tuple[str, str, str]:
    tm_text = "" if item.tm_celsius is None else f"{item.tm_celsius:.2f}"
    return (item.id, str(item.length), tm_text)
