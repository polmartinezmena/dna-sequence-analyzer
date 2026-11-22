"""Optional PNG plots for window GC and GC skew.

matplotlib is not installed with the base package. Importing this
module does not import matplotlib. The functions below import it only
when a plot is actually requested, and they explain how to install the
``plot`` extra if it is missing.

The figures are a quick look at the table, not a publication graphic.
"""

from collections import OrderedDict
from pathlib import Path

from dna_sequence_analyzer.windows import GCWindow, SkewWindow


def lab_template(mode: str) -> str:
    """Register and return a transparent Plotly theme for the current palette."""
    from plotly import graph_objects as go, io as pio

    from dna_sequence_analyzer.ui_theme import ACCENT, SERIES_2, STATUS, SURFACES

    name = "lab" if mode == "light" else "lab_dark"
    if name not in pio.templates:
        _, _, border, text, _ = SURFACES[mode]
        pio.templates[name] = go.layout.Template(layout={
            "font": {"family": "Inter, system-ui, sans-serif", "color": text},
            "xaxis": {"gridcolor": border},
            "yaxis": {"gridcolor": border},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)",
            "height": 320,
            "hovermode": "x unified",
            "legend": {"orientation": "h", "x": 0, "y": 1.12},
            "colorway": [ACCENT, SERIES_2, STATUS["Aviso"]],
            "separators": ",.",
        })
    return name


def write_gc_window_plot(rows: list[GCWindow], path: str | Path) -> None:
    """Plot GC% against the 1-based window start. One line per record."""
    plt = _pyplot()
    destination = _check_destination(path)
    fig, axis = plt.subplots(figsize=(8, 4))
    for record_id, group in _grouped(rows):
        points = [row for row in group if row.gc_percent is not None]
        axis.plot(
            [row.start for row in points],
            [row.gc_percent for row in points],
            marker="o",
            label=record_id,
        )
    axis.set_xlabel("Window start (1-based)")
    axis.set_ylabel("GC%")
    axis.set_title("GC content along the sequence")
    if rows:
        axis.legend()
    fig.tight_layout()
    fig.savefig(destination)
    plt.close(fig)


def write_skew_plot(rows: list[SkewWindow], path: str | Path) -> None:
    """Plot window skew and cumulative skew. One line per record."""
    plt = _pyplot()
    destination = _check_destination(path)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for record_id, group in _grouped(rows):
        defined = [row for row in group if row.gc_skew is not None]
        axes[0].plot(
            [row.start for row in defined],
            [row.gc_skew for row in defined],
            marker="o",
            label=record_id,
        )
        cumulative = [row for row in group if row.cumulative_skew is not None]
        axes[1].plot(
            [row.start for row in cumulative],
            [row.cumulative_skew for row in cumulative],
            marker="o",
            label=record_id,
        )
    axes[0].set_ylabel("GC skew")
    axes[0].set_title("(G - C) / (G + C) in each window")
    axes[1].set_xlabel("Window start (1-based)")
    axes[1].set_ylabel("Cumulative skew")
    axes[1].set_title("Running sum, restarted for each record")
    if rows:
        axes[0].legend()
    fig.tight_layout()
    fig.savefig(destination)
    plt.close(fig)


def _grouped(rows):
    """Group rows by id, keeping the order of the file."""
    groups: OrderedDict[str, list] = OrderedDict()
    for row in rows:
        groups.setdefault(row.id, []).append(row)
    return groups.items()


def _check_destination(path: str | Path) -> Path:
    destination = Path(path)
    if not destination.parent.exists():
        raise ValueError(f"Directory does not exist: {destination.parent}")
    return destination


def _pyplot():
    """Import pyplot, or raise ImportError with install instructions."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Plotting needs matplotlib. Install the optional extra with: "
            'pip install -e ".[plot]"'
        ) from exc
    return plt
