"""Shared presentation tokens for the workspace and its reports."""

BASE_COLORS = {
    "A": "#2E9E5B", "C": "#2F7FD1", "G": "#D08A1C",
    "T": "#D0453B", "U": "#D0453B", "other": "#8A94A6",
}
ACCENT = "#4C5BD4"
SERIES_2 = "#0E9AA7"
STATUS = {"OK": "#1F8A4C", "Aviso": "#B26B00", "Error": "#C2352B"}
SURFACES = {
    "light": ("#FFFFFF", "#F6F7F9", "#E2E6EC", "#141A23", "#5B6574"),
    "dark": ("#0F1318", "#161B22", "#2A313C", "#E6EAF0", "#9AA4B2"),
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


def stylesheet(theme=None):
    """Build the same stylesheet for Streamlit and standalone HTML reports."""
    names = ("surface", "panel", "border", "text", "muted")

    def variables(mode):
        return ";".join(f"--lab-{name}:{value}" for name, value in zip(names, SURFACES[mode]))

    palettes = f":root {{{variables(theme or 'light')}}}"
    if theme is None:
        palettes += f"@media(prefers-color-scheme:dark){{:root{{{variables('dark')}}}}}"
    return f"""<style>
{palettes}
:root {{
  --lab-accent:{ACCENT};--lab-series-2:{SERIES_2};
  --lab-ok:{STATUS['OK']};--lab-warning:{STATUS['Aviso']};--lab-error:{STATUS['Error']};
  --lab-space-1:4px;--lab-space-2:8px;--lab-space-3:12px;
  --lab-space-4:16px;--lab-space-5:24px;--lab-space-6:32px;
  --lab-radius:6px;--lab-rule:1px solid var(--lab-border);
  --lab-mono:"IBM Plex Mono",Consolas,monospace;
  color-scheme:{theme or 'light dark'};
}}
body {{background:var(--lab-surface);color:var(--lab-text);font:14px/1.5 Inter,system-ui,sans-serif}}
.block-container {{max-width:92rem;padding-top:3.2rem}}
.lab-record {{padding:var(--lab-space-4);border:var(--lab-rule);border-radius:var(--lab-radius);
  background:var(--lab-panel);overflow-wrap:anywhere}}
.lab-record h1 {{font:600 1.45rem/1.3 var(--lab-mono);padding:0;margin:0}}
.lab-record p {{margin:var(--lab-space-1) 0 0;color:var(--lab-muted)}}
.lab-run {{font: .75rem/1.5 var(--lab-mono);overflow-wrap:anywhere}}
.sequence-scroll {{max-width:100%;overflow-x:auto;font:14px/1.5 var(--lab-mono);margin-bottom:12px}}
.sequence-row {{display:grid;grid-template-columns:5rem repeat(60,1ch);min-width:max-content;gap:0 2px}}
.sequence-row span {{text-align:center;min-width:0}}
.sequence-row .track {{text-align:left;color:var(--lab-muted)}}
.sequence-row .tick {{text-align:left;min-width:max-content;position:relative;z-index:1}}
.sequence-row .base {{border-bottom:2px solid var(--base);color:var(--lab-text)}}
.sequence-row .complement {{opacity:.65}}
.sequence-row .met {{color:var(--lab-ok);font-weight:bold}}
.sequence-row .stop {{color:var(--lab-error);font-weight:bold}}
.sequence-row .continuation {{color:var(--lab-muted)}}
.lab-protein {{font-family:var(--lab-mono);white-space:pre-wrap}}
.lab-report {{max-width:80rem;margin:auto;padding:var(--lab-space-5);overflow-x:auto}}
.lab-report table {{width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums}}
.lab-report th,.lab-report td {{padding:var(--lab-space-2);border:var(--lab-rule);text-align:left}}
.lab-report th {{background:var(--lab-panel)}}
</style>"""
