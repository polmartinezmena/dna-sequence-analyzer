"""Command-line interface for the DNA sequence analyzer.

``gc`` prints length, GC content, ambiguous-base count and a flag for
extreme GC. ``gc-window`` and ``skew`` slide a window along each record.
``tm`` prints a nearest-neighbor melting temperature for short oligos.
Later stages will add more commands.
"""

import csv
import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from typing import TextIO

import typer

from dna_sequence_analyzer import __version__
from dna_sequence_analyzer.blast import BlastOptions, search
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
        "Analiza secuencias de ADN y ARN. "
        "Calcula GC y Tm, traduce, busca ORFs y consulta BLAST."
    ),
)


@app.command("orfs")
def orfs(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Archivo FASTA de nucleótidos."),
    min_aa: int = typer.Option(30, "--min-aa", min=1, help="Longitud mínima en aminoácidos."),
    table: int = typer.Option(1, "--table", help="Tabla del código genético de NCBI."),
    nested: bool = typer.Option(False, "--nested", help="Incluye inicios anidados en el mismo marco."),
    allow_no_stop: bool = typer.Option(False, "--allow-no-stop", help="Incluye ORFs abiertos al final."),
    alt_starts: bool = typer.Option(
        False, "--alt-starts", help="Acepta también los codones de inicio alternativos de la tabla NCBI elegida"
    ),
    output_format: str = typer.Option("csv", "--format", help="fasta, csv o gff."),
) -> None:
    """Busca ORFs en las dos hebras y exporta coordenadas inclusivas desde 1."""
    try:
        hits = [(record.id, hit) for record in read_fasta(fasta_file) for hit in find_orfs(
            str(record.seq), min_aa=min_aa, table=table, include_nested=nested,
            allow_no_stop=allow_no_stop, alternative_starts=alt_starts,
        )]
        typer.echo(export_orfs(hits, output_format), nl=False)
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc


@app.command("blast")
def blast(
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Un registro FASTA público."),
    program: str = typer.Option("blastn", "--program", help="blastn para nucleótidos; blastp para proteínas."),
    db: str | None = typer.Option(None, "--db", help="nt, nr o refseq_rna según el programa."),
    evalue: float = typer.Option(0.01, "--evalue", help="Valor E máximo."),
    max_hits: int = typer.Option(20, "--max-hits", min=1, help="Número máximo de coincidencias."),
    megablast: bool = typer.Option(True, "--megablast/--no-megablast", help="Activa MegaBLAST para blastn."),
    timeout: float = typer.Option(300, "--timeout", help="Tiempo máximo en segundos."),
) -> None:
    """Envía un FASTA público a NCBI y escribe los resultados como CSV."""
    try:
        records = read_fasta(fasta_file)
        if len(records) != 1:
            raise ValueError("La consulta BLAST requiere un único registro FASTA.")
        record = records[0]
        options = BlastOptions(program, db or ("nr" if program == "blastp" else "nt"), evalue, max_hits, megablast)
        hits = search(str(record.seq), options, email=os.getenv("NCBI_EMAIL"), timeout=timeout)
    except (OSError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(("accession", "description", "identity_percent", "coverage_percent", "evalue", "bitscore"))
    writer.writerows((hit.accession, hit.description, hit.identity, hit.coverage, hit.evalue,
                      hit.bitscore) for hit in hits)


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
        help="Muestra la versión del programa.",
    ),
) -> None:
    """Analiza secuencias de ADN desde la línea de comandos."""


@app.command("gc")
def gc(
    fasta_file: Path = typer.Argument(
        ...,
        metavar="FASTA",
        help="Archivo FASTA con uno o más registros de ADN o ARN.",
    ),
    as_csv: bool = typer.Option(
        False,
        "--csv",
        help="Emite la tabla como CSV.",
    ),
) -> None:
    """Muestra longitud, porcentaje de GC y bases ambiguas por registro.

    El cálculo cuenta G, C y S como GC. A, T, U y W cuentan como bases
    conocidas que no son GC. Los extremos por debajo de 25 % y por encima
    de 65 % aparecen como avisos."""
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
        help="Archivo FASTA de oligos. Omítelo si usas --sequence.",
    ),
    sequence: str | None = typer.Option(
        None,
        "--sequence",
        "-s",
        help="Un oligo con A, C, G, T o U, introducido directamente.",
    ),
    na: float = typer.Option(50, "--na", help="Concentración de Na+ (mM)."),
    k: float = typer.Option(0, "--k", help="Concentración de K+ (mM)."),
    tris: float = typer.Option(0, "--tris", help="Concentración de Tris (mM)."),
    mg: float = typer.Option(0, "--mg", help="Concentración de Mg2+ (mM)."),
    dntps: float = typer.Option(0, "--dntps", help="Concentración de dNTP (mM)."),
    dnac1: float = typer.Option(
        25,
        "--dnac1",
        help="Concentración de la hebra más abundante (nM).",
    ),
    dnac2: float = typer.Option(
        25,
        "--dnac2",
        help="Concentración de la hebra menos abundante (nM).",
    ),
    saltcorr: int = typer.Option(
        5,
        "--saltcorr",
        help="Corrección de sal (0–7). 5 usa SantaLucia 1998. 0 la desactiva.",
        min=0,
        max=7,
    ),
    self_complementary: bool = typer.Option(
        False,
        "--self-complementary",
        help="Trata el oligo como autocomplementario e ignora dnac2.",
    ),
    as_csv: bool = typer.Option(
        False,
        "--csv",
        help="Emite id,length,tm_celsius como CSV.",
    ),
) -> None:
    """Estima la temperatura de fusión de oligos cortos por vecinos más próximos.

    Pasa un FASTA o usa --sequence para una secuencia individual."""
    if fasta_file is None and sequence is None:
        typer.echo("Error: Indica un archivo FASTA o --sequence.", err=True)
        raise typer.Exit(code=1)
    if fasta_file is not None and sequence is not None:
        typer.echo("Error: Indica un archivo FASTA o --sequence como única entrada.", err=True)
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
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Archivo FASTA."),
    window: int = typer.Option(..., "--window", min=1, help="Longitud de la ventana en bases."),
    step: int = typer.Option(
        ...,
        "--step",
        min=1,
        help="Bases entre inicios de ventana. Igual a --window evita solapamiento.",
    ),
    as_csv: bool = typer.Option(False, "--csv", help="Emite la tabla como CSV."),
    plot: Path | None = typer.Option(
        None,
        "--plot",
        help="Escribe el perfil de GC en PNG. Requiere el extra plot.",
    ),
) -> None:
    """Muestra GC y CpG observado/esperado para cada ventana deslizante.

    Las coordenadas empiezan en 1 e incluyen ambos extremos."""
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
    fasta_file: Path = typer.Argument(..., metavar="FASTA", help="Archivo FASTA."),
    window: int = typer.Option(..., "--window", min=1, help="Longitud de la ventana en bases."),
    step: int = typer.Option(..., "--step", min=1, help="Bases entre inicios de ventana."),
    as_csv: bool = typer.Option(False, "--csv", help="Emite la tabla como CSV."),
    plot: Path | None = typer.Option(
        None,
        "--plot",
        help="Escribe el sesgo local y acumulado en PNG. Requiere el extra plot.",
    ),
) -> None:
    """Muestra el sesgo GC por ventana y su suma acumulada."""
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
    port: int = typer.Option(8501, "--port", min=1, help="Puerto de la app web local."),
) -> None:
    """Abre la app local en el navegador. Requiere el extra web."""
    if importlib.util.find_spec("streamlit") is None:
        typer.echo(
            'Error: La app web requiere Streamlit. Instálalo con pip install -e ".[web]"',
            err=True,
        )
        raise typer.Exit(code=1)
    entry = streamlit_entrypoint()
    if not entry.is_file():
        typer.echo(f"Error: No se encontró la app web en {entry}.", err=True)
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
    from dna_sequence_analyzer.ui_theme import ACCENT, APP_BG, BORDER, CARD, TEXT

    return [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(streamlit_entrypoint()),
        "--theme.base", "light",
        "--theme.primaryColor", ACCENT,
        "--theme.backgroundColor", APP_BG,
        "--theme.secondaryBackgroundColor", CARD,
        "--theme.textColor", TEXT,
        "--theme.borderColor", BORDER,
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
    typer.echo(f"Archivo guardado en {path}", err=True)


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
