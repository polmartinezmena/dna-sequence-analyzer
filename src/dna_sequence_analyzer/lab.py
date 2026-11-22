"""Quality checks and reproducible exports for the web interface."""
from datetime import datetime, timezone
from hashlib import sha256
from uuid import uuid4
import csv
import io
import json

from jinja2 import Environment

from . import __version__
from .ui_theme import format_es


def quality(records, stats, error=None, *, comments=0):
    if error:
        return [("Error", error)]
    checks = []
    for record, item in zip(records, stats):
        bases = str(record.seq).upper()
        invalid = sorted(set(bases) - set("ACGTURYSWKMBDHVN"))
        if not bases:
            checks.append(("Error", f"{record.id} · Secuencia vacía."))
        if invalid:
            checks.append(("Error", f"{record.id} · Caracteres inválidos: {', '.join(invalid)}."))
        ambiguous = sum(base in "RYSWKMBDHVN" for base in bases)
        if ambiguous:
            checks.append(("Aviso", f"{record.id} · {ambiguous} bases ambiguas."))
        if item.gc_flag:
            checks.append(("Aviso", f"{record.id} · GC extremo ({item.gc_percent:.1f} %)."))
    if comments:
        checks.append(("Aviso", f"se ignoraron {format_es(comments)} l?neas de comentario antes del primer registro"))
    return checks or [("OK", "Muestra sin incidencias.")]


def run_metadata(records, parameters, previous=None):
    sample = [(r.id, str(r.seq)) for r in records]
    digest = sha256(json.dumps(sample, ensure_ascii=True).encode()).hexdigest()
    values = {"SHA-256 de la muestra": digest, "Versión": __version__, **parameters}
    if previous and all(previous.get(k) == v for k, v in values.items()) and len(previous) == len(values) + 2:
        return previous
    return {"ID": uuid4().hex, "Fecha (UTC)": datetime.now(timezone.utc).isoformat(), **values}


def result_csv(rows):
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def report(metadata, checks, results):
    template = Environment(autoescape=True).from_string("""<!doctype html>
<html lang="es"><meta charset="utf-8"><title>Informe de análisis de ADN</title>
<h1>Informe de análisis de ADN</h1><h2>Metadatos del análisis</h2><dl>
{% for k,v in metadata.items() %}<dt><b>{{ k }}</b></dt><dd>{{ v }}</dd>{% endfor %}</dl>
<h2>Control de calidad</h2><ul>{% for level,msg in checks %}<li>{{ level }} · {{ msg }}</li>{% endfor %}</ul>
{% for title,rows in results.items() %}<h2>{{ title }}</h2>{% if rows %}<table border="1">
<tr>{% for k in rows[0] %}<th>{{ k }}</th>{% endfor %}</tr>
{% for row in rows %}<tr>{% for v in row.values() %}<td>{{ v if v is not none else "—" }}</td>{% endfor %}</tr>{% endfor %}
</table>{% endif %}{% endfor %}</html>""")
    return template.render(metadata=metadata, checks=checks, results=results)
