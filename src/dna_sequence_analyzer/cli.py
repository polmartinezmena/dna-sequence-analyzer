"""Command-line interface for the DNA sequence analyzer.

``gc`` prints length, GC content, ambiguous-base count and a flag for
extreme GC. ``gc-window`` and ``skew`` slide a window along each record.
``tm`` prints a nearest-neighbor melting temperature for short oligos.
Later stages will add more commands.
"""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import TextIO

import typer

from dna_sequence_analyzer import __version__
from dna_sequence_analyzer.gc import GC_HIGH_PERCENT, GC_LOW_PERCENT, analyze_fasta, read_fasta
from dna_sequence_analyzer.orfs import export_orfs, find_orfs
from dna_sequence_analyzer.translate import translate_frames
from dna_sequence_analyzer.plots import write_gc_window_plot, write_skew_plot
from dna_sequence_analyzer.tm import (
    OLIGO_GUIDANCE_LENGTH,
    MeltingConditions,
    MeltingRecord,
    analyze_melting,
    melting_record_for_sequence,
)
from dna_sequence_analyzer.windows import gc_windows, skew_windows

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "DNA sequence analysis toolkit. "
        "Reads FASTA and reports GC content, window composition, GC skew, and oligo Tm."
    ),
)


@app.command("orfs")
def orfs(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Archivo FASTA de nucleótidos."),
    min_aa: int = typer.Option(30, "--min-aa", min=1, help="Longitud mínima en aminoácidos."),
    table: int = typer.Option(1, "--table", help="Tabla del código genético de NCBI."),
    nested: bool = typer.Option(False, "--nested", help="Incluye inicios anidados en el mismo marco."),
    allow_no_stop: bool = typer.Option(False, "--allow-no-stop", help="Incluye ORFs abiertos al final."),
    output_format: str = typer.Option("csv", "--format", help="fasta, csv o gff."),
) -> None:
    """Busca ORFs en las dos hebras y exporta coordenadas inclusivas desde 1."""
    try:
        hits = [(record.id, hit) for record in read_fasta(fasta_file) for hit in find_orfs(
            str(record.seq), min_aa=min_aa, table=table, include_nested=nested,
            allow_no_stop=allow_no_stop,
        )]
        typer.echo(export_orfs(hits, output_format), nl=False)
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc


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
        help="Print the table as CSV instead of aligned text.",
    ),
) -> None:
    """Print length, GC content and ambiguous-base count for each record.

    GC% counts G, C and S as GC, and A, T, U and W as known bases that
    are not GC. Other ambiguity codes, such as N, are counted in the
    ambiguous column and excluded from the percentage. U means the
    record is reported as RNA. GC below 25% or above 65% is flagged,
    because those extremes often bias PCR and Illumina sequencing.
    """
    try:
        stats = analyze_fasta(fasta_file)
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    header = ("id", "length", "gc_percent", "ambiguous", "type", "flag")
    rows = [
        (
            item.id,
            str(item.length),
            "" if item.gc_percent is None else f"{item.gc_percent:.2f}",
            str(item.ambiguous_bases),
            item.nucleic_acid,
            item.gc_flag or "",
        )
        for item in stats
    ]
    _emit_table(header, rows, as_csv, na_columns=(2,))
    if any(item.gc_percent is None for item in stats):
        typer.echo(
            "Note: GC% is NA when a record has no bases that can be counted "
            "(A, C, G, T, U, S or W).",
            err=True,
        )
    if any(item.nucleic_acid in {"rna", "mixed"} for item in stats):
        typer.echo(
            "Note: U marks RNA and counts as a non-GC base, like T. "
            "RNA secondary structure and G-U wobble pairs are not predicted.",
            err=True,
        )
    if any(item.gc_flag for item in stats):
        typer.echo(
            f"Note: GC below {GC_LOW_PERCENT:.0f}% is flagged low-gc and GC above "
            f"{GC_HIGH_PERCENT:.0f}% is flagged high-gc. Extreme GC can bias "
            "PCR amplification and Illumina sequencing.",
            err=True,
        )


@app.command("translate")
def translate(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Archivo FASTA de entrada."),
    table: int = typer.Option(1, "--table", help="Número de tabla del código genético NCBI."),
    frame: str = typer.Option("all", "--frame", help="all, +1, +2, +3, -1, -2 o -3."),
) -> None:
    """Resume la traducción de cada secuencia por marco de lectura."""
    rows = []
    try:
        for record in read_fasta(fasta_file):
            for item in translate_frames(str(record.seq), table, frame):
                rows.append((record.id, item.frame, str(len(item.protein) - item.stops),
                             str(item.stops), str(item.longest_stretch)))
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(_format_table(("id", "marco", "longitud_aa", "paradas", "tramo_máximo"), rows))


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
    if any(item.contains_uracil for item in records):
        typer.echo(
            "Note: U was read as T so the DNA/DNA nearest-neighbor table can score it. "
            "This is not an RNA folding temperature, and G-U wobble pairs are not modeled.",
            err=True,
        )


@app.command("gc-window")
def gc_window(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="FASTA file."),
    window: int = typer.Option(..., "--window", min=1, help="Window length in bases."),
    step: int = typer.Option(
        ...,
        "--step",
        min=1,
        help="Bases between window starts. Equal to --window means no overlap.",
    ),
    as_csv: bool = typer.Option(False, "--csv", help="Print CSV instead of a text table."),
    plot: Path | None = typer.Option(
        None,
        "--plot",
        help="Write a PNG of GC% along the sequence. Needs the plot extra.",
    ),
) -> None:
    """Print GC% and CpG observed/expected for each sliding window.

    Coordinates are 1-based and include both ends. A window that would
    pass the end of a record is omitted. CpG O/E is
    (CG dinucleotides * length) / (C count * G count). It is one of the
    numbers used to notice CpG islands; this command does not call
    islands by itself.
    """
    try:
        rows = gc_windows(fasta_file, window, step)
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    if not rows:
        typer.echo(
            "Note: no complete window fitted. Every record is shorter than --window.",
            err=True,
        )
        raise typer.Exit(code=0)

    header = ("id", "start", "end", "gc_percent", "cpg_oe", "ambiguous")
    table = [
        (
            item.id,
            str(item.start),
            str(item.end),
            "" if item.gc_percent is None else f"{item.gc_percent:.2f}",
            "" if item.cpg_oe is None else f"{item.cpg_oe:.4f}",
            str(item.ambiguous_bases),
        )
        for item in rows
    ]
    _emit_table(header, table, as_csv, na_columns=(3, 4))
    _maybe_plot(plot, lambda destination: write_gc_window_plot(rows, destination))


@app.command("skew")
def skew(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="FASTA file."),
    window: int = typer.Option(..., "--window", min=1, help="Window length in bases."),
    step: int = typer.Option(..., "--step", min=1, help="Bases between window starts."),
    as_csv: bool = typer.Option(False, "--csv", help="Print CSV instead of a text table."),
    plot: Path | None = typer.Option(
        None,
        "--plot",
        help="Write a PNG of window skew and cumulative skew. Needs the plot extra.",
    ),
) -> None:
    """Print GC skew and cumulative skew for each sliding window.

    Skew is (G - C) / (G + C). Cumulative skew is the running sum of
    those values and restarts on every record. On many bacterial
    chromosomes the cumulative minimum sits near the replication origin
    and the maximum near the terminus. That is a clue, not a proof.
    """
    try:
        rows = skew_windows(fasta_file, window, step)
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    if not rows:
        typer.echo(
            "Note: no complete window fitted. Every record is shorter than --window.",
            err=True,
        )
        raise typer.Exit(code=0)

    header = ("id", "start", "end", "gc_skew", "cumulative_skew")
    table = [
        (
            item.id,
            str(item.start),
            str(item.end),
            "" if item.gc_skew is None else f"{item.gc_skew:.4f}",
            "" if item.cumulative_skew is None else f"{item.cumulative_skew:.4f}",
        )
        for item in rows
    ]
    _emit_table(header, table, as_csv, na_columns=(3, 4))
    _maybe_plot(plot, lambda destination: write_skew_plot(rows, destination))


@app.command("web")
def web(
    port: int = typer.Option(8501, "--port", min=1, help="Port for the local web app."),
) -> None:
    """Open the browser app (Spanish interface).

    Install it first with: pip install -e ".[web]"
    """
    if importlib.util.find_spec("streamlit") is None:
        typer.echo(
            'Error: the web app needs Streamlit. Install it with: pip install -e ".[web]"',
            err=True,
        )
        raise typer.Exit(code=1)
    entry = streamlit_entrypoint()
    if not entry.is_file():
        typer.echo(f"Error: cannot find the web app at {entry}.", err=True)
        raise typer.Exit(code=1)
    raise typer.Exit(subprocess.call(web_launch_command(port)))


def streamlit_entrypoint() -> Path:
    """File that ``streamlit run`` should open.

    Prefer ``app/streamlit_app.py`` in a repository checkout, because
    that is also the file Streamlit Community Cloud runs. After a
    normal install, fall back to the module shipped inside the package.
    """
    repo_launcher = Path(__file__).resolve().parents[2] / "app" / "streamlit_app.py"
    if repo_launcher.is_file():
        return repo_launcher
    return Path(__file__).resolve().parent / "webapp.py"


def web_launch_command(port: int = 8501) -> list[str]:
    """Command list used by ``dna-analyzer web``. Exposed for tests."""
    return [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(streamlit_entrypoint()),
        "--server.port",
        str(port),
    ]


def main() -> None:
    """Entry point used by the ``dna-analyzer`` console script."""
    app()


def _format_optional(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.2f}"


def _emit_table(
    header: tuple[str, ...],
    rows: list[tuple[str, ...]],
    as_csv: bool,
    na_columns: tuple[int, ...] = (),
) -> None:
    """Print rows. CSV keeps undefined numbers blank. Text shows NA."""
    if as_csv:
        _write_csv(header, rows, sys.stdout)
        return
    display = []
    for row in rows:
        cells = list(row)
        for index in na_columns:
            if cells[index] == "":
                cells[index] = "NA"
        display.append(tuple(cells))
    typer.echo(_format_table(header, display))


def _maybe_plot(path: Path | None, writer) -> None:
    if path is None:
        return
    try:
        writer(path)
    except (ImportError, ValueError, OSError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"Wrote {path}", err=True)


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
    # An empty last cell is padding. Drop that trailing space so a row
    # with no flag does not end in invisible characters.
    return "\n".join(line.rstrip() for line in lines)


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
