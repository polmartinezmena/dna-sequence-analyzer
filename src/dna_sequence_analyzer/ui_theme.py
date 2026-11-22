"""Presentation tokens, formatting and the Streamlit workspace stylesheet."""

import re

ACCENT = "#E08A2E"
ACCENT_TEXT = "#91480E"
INDIGO_SOFT = "#FAECD9"
APP_BG = "#F6F4EF"
CARD = "#FFFFFF"
BORDER = "#ECE8E1"
TEXT = "#28231E"
MUTED = "#6B6258"
DARK_PANEL = "#1F1B16"
DARK_GRID = "#393129"
DARK_TEXT = "#F0E9DF"
DARK_MUTED = "#C1B4A6"
CHART_LINE = "#22B8CF"
REFERENCE = "#ACA094"
SERIES_2 = "#49837A"
OK = "#23784B"
AVISO = "#9B570F"
ERROR = "#C2352B"
STATUS = {"OK": OK, "Aviso": AVISO, "Error": ERROR}
BASE_COLORS = {
    "A": "#EAA762", "C": "#83CAD4", "G": "#8BC494",
    "T": "#EAD176", "U": "#EAD176", "other": "#716C68",
}


def format_es(value, decimals=0, unit=""):
    """Format a number for display, keeping undefined values distinct from zero."""
    if value is None:
        return "—"
    number = f"{value:,.{decimals}f}".translate(str.maketrans(",.", ".,"))
    return f"{number} {unit}" if unit else number


def spanish_table(rows, decimals=None):
    """Format displayed values while preserving numeric types for sorting."""
    import pandas as pd

    styler = pd.DataFrame(rows).style.format(decimal=",", thousands=".", precision=0, na_rep="—")
    for column, precision in (decimals or {}).items():
        styler.format(subset=[column], decimal=",", thousands=".", precision=precision, na_rep="—")
    return styler


def stylesheet(*, active_view=None):
    """Build a light stylesheet for the application and its standalone reports."""
    tokens = {
        "bg": APP_BG, "card": CARD, "border": BORDER, "text": TEXT, "muted": MUTED,
        "accent": ACCENT, "accent-text": ACCENT_TEXT, "active": INDIGO_SOFT,
        "dark-panel": DARK_PANEL, "dark-grid": DARK_GRID, "dark-text": DARK_TEXT,
        "dark-muted": DARK_MUTED, "ok": OK, "warning": AVISO, "error": ERROR,
    }
    palette = ";".join(f"--lab-{name}:{color}" for name, color in tokens.items())
    active = ""
    if active_view:
        key = re.sub(r"[^a-zA-Z0-9_-]", "-", active_view)
        active = (
            f'.st-key-nav-{key} button {{background:var(--lab-card);color:var(--lab-accent-text);'
            'box-shadow:inset 3px 0 var(--lab-accent);font-weight:650}'
        )
    # Fix the number of cells per row at each available width so CSS can reveal the row-start tick.
    tile_wrap = "\n".join(
        f"@container base-ruler (min-width:{20 * cols}px){maximum}{{"
        f".tile-grid{{width:{20 * cols}px}}"
        f".tile-cell:nth-child({cols}n+1) .tile-tick{{visibility:visible}}}}"
        for cols in range(2, 61)
        for maximum in (f" and (max-width:{20 * (cols + 1) - 0.01:.2f}px)" if cols < 60 else "",)
    )
    return f"""<style>
:root {{{palette};--lab-mono:"IBM Plex Mono",Consolas,monospace;
  --lab-radius:18px;color-scheme:light}}
{_CSS}
{active}
{tile_wrap}
</style>"""


_CSS = """
body {background:var(--lab-bg);color:var(--lab-text);font:14px/1.5 Inter,system-ui,sans-serif}
::selection {background:var(--lab-active);color:var(--lab-text)}
[data-testid="stAppViewContainer"] {background:var(--lab-bg)}
header[data-testid="stHeader"] {background:transparent;height:2.7rem}
[data-testid="stSidebar"] {background:var(--lab-bg);border:none}
[data-testid="stSidebarContent"] {border:none}
[data-testid="stSidebarUserContent"] {padding:2.8rem 1rem 1.1rem}
.block-container {max-width:104rem;padding:2.9rem clamp(1rem,2.5vw,2.6rem) 3rem}
.topbar {height:56px;display:flex;align-items:center;gap:12px;padding:0 20px;
  background:var(--lab-card);border-bottom:1px solid var(--lab-border);
  border-radius:var(--lab-radius) var(--lab-radius) 0 0;color:var(--lab-muted)}
.topbar strong {color:var(--lab-text);font-weight:700}
.topbar-separator {color:#AA9F92}
[data-testid="stHeading"] h1 {font-size:1.6rem;letter-spacing:-.025em;font-weight:700;padding:.4rem 0}
.section-label {font-size:.8rem;font-weight:750;line-height:1.4;letter-spacing:.09em;
  text-transform:uppercase;color:var(--lab-accent-text);margin:0;padding:0}
.st-key-muestra, .st-key-navegacion {gap:.25rem}
.sample-tree {display:flex;align-items:baseline;gap:.45rem;margin:.35rem 0 .1rem;overflow-wrap:anywhere}
.sample-tree .file-icon {font:1.1rem "Material Symbols Rounded";color:var(--lab-muted)}
.st-key-muestra [data-testid="stButton"] button, .st-key-navegacion button {justify-content:flex-start}
.st-key-muestra [class*="st-key-registro-"] {padding-left:.9rem}
.st-key-muestra [class*="st-key-registro-"] button {min-height:2.15rem}
.st-key-muestra [class*="st-key-registro-"] code {font:500 .88rem var(--lab-mono);
  padding:0;background:transparent;color:inherit}
.st-key-muestra :is(.st-key-abrir-muestra,.st-key-pegar-muestra) button {padding:.3rem .55rem;min-height:2rem}
.st-key-navegacion button {min-height:2.55rem;color:var(--lab-muted);padding-left:1rem}
.st-key-navegacion button:hover {background:var(--lab-active);color:var(--lab-text)}
.st-key-navegacion :is([class*="Oligo"],.st-key-nav-Aprende) button p {font-size:.88rem}
.st-key-navegacion hr {margin:.65rem 0}
.lab-run {font:.75rem/1.55 var(--lab-mono);overflow-wrap:anywhere;color:var(--lab-muted);padding:.7rem .2rem}
:is(.st-key-cabecera-muestra,.st-key-panel-lateral,.st-key-parametros,[class*="st-key-tarjeta-"]){
  background:var(--lab-card);border:1px solid var(--lab-border);border-radius:var(--lab-radius);
  padding:24px;min-width:0}
.sample-identity {display:flex;align-items:center;gap:12px;min-width:0}
.sample-avatar {display:grid;place-items:center;flex:0 0 42px;height:42px;
  background:var(--lab-active);border-radius:50%;font-weight:750;color:var(--lab-accent-text)}
.sample-date {font-size:.72rem;letter-spacing:.08em;text-transform:uppercase;color:var(--lab-muted);margin:0}
.sample-filename {font-size:1.45rem;line-height:1.15;letter-spacing:-.025em;
  font-weight:750;overflow-wrap:anywhere;margin:.25rem 0 0}
.sample-stat {border-left:1px solid var(--lab-border);padding-left:1rem}
.sample-stat strong {display:block;font-size:1.35rem;font-weight:750;line-height:1.2;
  white-space:nowrap;font-variant-numeric:tabular-nums}
.sample-stat small {display:block;color:var(--lab-muted);margin-top:.2rem}
.st-key-cabecera-muestra [data-testid="stColumn"]:nth-child(4) strong {color:var(--lab-accent-text)}
button[data-testid^="stBaseButton-"] {border-radius:999px}
button[data-testid^="stBaseButton-primary"]:not(:disabled) {background:var(--lab-accent);
  color:var(--lab-text);border-color:var(--lab-accent);font-weight:650}
button[data-testid^="stBaseButton-primary"]:hover:not(:disabled) {background:#F0A657;color:var(--lab-text)}
button[data-testid^="stBaseButton-secondary"]:not(:disabled) {background:var(--lab-card);
  color:var(--lab-text);border:1px solid var(--lab-border)}
button[data-testid^="stBaseButton-secondary"]:hover:not(:disabled) {border-color:var(--lab-accent);
  color:var(--lab-accent-text)}
.st-key-cabecera-muestra button {white-space:nowrap}
.st-key-cabecera-muestra button p {font-size:.88rem}
.st-key-cabecera-muestra .st-key-importar-cabecera button {color:var(--lab-accent-text)}
.st-key-panel-oscuro {background:var(--lab-dark-panel);border-radius:12px;padding:16px 20px;
  min-width:0;color:var(--lab-dark-text);--lab-text:var(--lab-dark-text);--lab-muted:var(--lab-dark-muted)}
.st-key-panel-oscuro :is([data-testid="stWidgetLabel"] p,[data-testid="stMarkdownContainer"] p,
  [data-testid="stAlert"] p) {color:var(--lab-dark-text)}
.st-key-panel-oscuro [data-testid="stAlert"] {background:var(--lab-dark-grid)}
.st-key-panel-oscuro [role="tab"] p {color:var(--lab-dark-muted)}
.st-key-panel-oscuro [role="tab"][aria-selected="true"] p {color:#FFF}
.st-key-panel-oscuro [data-baseweb="tab-highlight"] {background:#FFF}
.st-key-panel-oscuro [data-baseweb="tab-border"] {background:var(--lab-dark-grid)}
.st-key-panel-oscuro [data-baseweb="input"],
.st-key-panel-oscuro [data-baseweb="base-input"] {background:var(--lab-dark-grid);border-color:var(--lab-dark-grid)}
.st-key-panel-oscuro .stNumberInput :is(input,button) {color:var(--lab-dark-text);background:transparent}
.base-tiles {container:base-ruler / inline-size;min-width:0;margin:.7rem 0 1rem}
.tile-ruler-label {color:var(--lab-dark-muted);font-size:.82rem;margin:0 0 .5rem}
.tile-record {display:inline-block;margin:0 0 .55rem;padding:2px 8px;border-radius:4px;
  font:.86rem/1.6 var(--lab-mono);background:var(--lab-dark-grid);color:var(--lab-dark-text)}
.tile-grid {display:flex;flex-wrap:wrap;column-gap:2px;row-gap:.5rem;max-width:100%}
.tile-cell {display:grid;grid-template-rows:16px 22px;width:18px;flex:0 0 18px}
.tile-tick {visibility:hidden;min-width:max-content;white-space:nowrap;color:var(--lab-dark-muted);
  font:10px/16px var(--lab-mono);position:relative;z-index:1}
.tile-tick.mark {visibility:visible}
.tile-base {display:grid;place-items:center;background:var(--base);border-radius:4px;
  color:#1F1B16;font:650 12px var(--lab-mono)}
.tile-ambiguous {color:#FFF}
.orf-tracks {display:grid;gap:6px;padding:.65rem 0}
.orf-track {display:grid;grid-template-columns:2.5rem 1fr;align-items:center;
  color:var(--lab-dark-text);font:.85rem var(--lab-mono)}
.orf-track-rail {position:relative;height:11px;background:var(--lab-dark-grid);border-radius:4px}
.orf-track-bar {position:absolute;inset-block:0;background:var(--lab-accent);border-radius:4px}
.qc-list {list-style:none;display:grid;gap:.45rem;margin:.7rem 0 1.1rem;padding:0}
.qc-list li {display:flex;gap:.55rem;align-items:baseline;color:var(--lab-text);font-size:.9rem;
  line-height:1.4;overflow-wrap:anywhere}
.qc-icon {font-family:"Material Symbols Rounded";font-feature-settings:"liga";font-size:1.08rem;color:var(--lab-ok)}
.qc-aviso .qc-icon {color:var(--lab-warning)}
.qc-error .qc-icon {color:var(--lab-error)}
.record-heading {border-top:1px solid var(--lab-border);padding-top:1rem;margin-top:.8rem}
.record-details {margin:.7rem 0 0;font-size:.9rem}
.record-details div {display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:.35rem;padding:.2rem 0}
.record-details dt {color:var(--lab-muted)}
.record-details dd {margin:0;text-align:right;font-weight:600;overflow-wrap:anywhere}
.record-id {font-family:var(--lab-mono)}
.region-dash {text-align:center;color:var(--lab-muted);margin-top:1.2rem}
.st-key-form-region [data-testid="stNumberInputStepDown"], .st-key-form-region [data-testid="stNumberInputStepUp"] {
  display:none}
.sequence-scroll {max-width:100%;overflow-x:auto;font:14px/1.5 var(--lab-mono);margin-bottom:12px}
.sequence-row {display:grid;grid-template-columns:5rem repeat(60,1ch);min-width:max-content;gap:0 2px}
.sequence-row span {text-align:center;min-width:0}
.sequence-row .track {text-align:left;color:var(--lab-muted)}
.sequence-row .tick {text-align:left;min-width:max-content;position:relative;z-index:1}
.sequence-row .base {border-bottom:2px solid var(--base);color:var(--lab-text)}
.sequence-row .complement {opacity:.75}
.sequence-row .met {color:var(--lab-ok);font-weight:bold}
.sequence-row .stop {color:var(--lab-error);font-weight:bold}
.st-key-panel-oscuro .sequence-row .met {color:#89D9A9}
.st-key-panel-oscuro .sequence-row .stop {color:#FFA99E}
.lab-protein {font-family:var(--lab-mono);white-space:pre-wrap}
.lab-report {max-width:80rem;margin:auto;padding:24px;overflow-x:auto}
.lab-report table {width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums}
.lab-report th,.lab-report td {padding:8px;border:1px solid var(--lab-border);text-align:left}
.lab-report th {background:var(--lab-card)}
@media(max-width:1240px) {
  .st-key-cabecera-muestra [data-testid="stColumn"] [data-testid="stHorizontalBlock"] {flex-wrap:wrap;gap:.3rem}
  .st-key-cabecera-muestra [data-testid="stColumn"] [data-testid="stColumn"] {min-width:max-content}
}
@media(max-width:850px) {
  :is(.st-key-cabecera-muestra,.st-key-panel-lateral,.st-key-parametros,[class*="st-key-tarjeta-"]){padding:16px}
  .st-key-cabecera-muestra [data-testid="stHorizontalBlock"] {flex-wrap:wrap;gap:.6rem}
  .st-key-cabecera-muestra [data-testid="stColumn"] {min-width:max-content}
}
@media(max-width:720px) {
  .st-key-espacio-vista > [data-testid="stHorizontalBlock"] {flex-wrap:wrap}
  .st-key-espacio-vista > [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {min-width:100%}
  .topbar {padding:0 12px;gap:8px}
}
"""
