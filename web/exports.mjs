import {orfs} from './core.mjs';

const html = value => String(value).replace(/[&<>"']/g, char => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
}[char]));
function csv(value) {
  const text = String(value);
  const safe = typeof value === 'string' && /^\s*[=+\-@\t\r\n]/.test(text) ? "'" + text : text;
  return '"' + safe.replaceAll('"', '""') + '"';
}
const gffId = value => encodeURIComponent(value || 'sequence').replace(/[!'()*]/g,
  char => '%' + char.charCodeAt(0).toString(16).toUpperCase());

export function exportSample(format, sample) {
  const {name, records, stats, table, orfOptions} = sample;
  const rows = stats.map(row => [row.id, row.length, row.gc_percent ?? '', row.ambiguous_bases,
    row.nucleic_acid, row.gc_flag ?? '']);
  if (format === 'csv') {
    const header = 'id,length,gc_percent,ambiguous_bases,type,gc_flag';
    return {type: 'text/csv;charset=utf-8', text: [header, ...rows.map(row => row.map(csv).join(','))].join('\n') + '\n'};
  }
  if (format === 'html') {
    const formatNumber = new Intl.NumberFormat('es-ES', {maximumFractionDigits: 2});
    const values = rows.map(row => '<tr>' + row.map(value => '<td>' +
      html(typeof value === 'number' ? formatNumber.format(value) : value) + '</td>').join('') + '</tr>');
    return {type: 'text/html;charset=utf-8', text: '<!doctype html><html lang="es"><meta charset="utf-8">' +
      '<title>Informe FASTA</title><style>body{font:16px system-ui;margin:2rem}th,td{padding:.45rem;text-align:left}' +
      'table{border-collapse:collapse}td{border-top:1px solid #ddd}</style><h1>' + html(name) + '</h1>' +
      '<h2>Registros</h2><table><thead><tr><th>Identificador</th><th>Longitud</th><th>GC (%)</th>' +
      '<th>Ambiguas</th><th>Tipo</th><th>Aviso</th></tr></thead><tbody>' + values.join('') + '</tbody></table></html>'};
  }
  if (format !== 'fasta' && format !== 'gff3') throw new Error('Formato desconocido.');
  const hits = records.flatMap(record => {
    try {
      return orfs(record.sequence, {...orfOptions, table}).map(hit => ({...hit, id: record.id}));
    } catch (error) {
      throw new Error((record.id || '(sin ID)') + ': ' + error.message);
    }
  });
  if (!hits.length && format === 'fasta') throw new Error('La muestra no tiene ORFs para exportar.');
  if (format === 'fasta') {
    return {type: 'text/plain;charset=utf-8', text: hits.map((hit, index) => '>orf' + (index + 1) +
      ' ' + hit.id + ' ' + hit.frame + ' ' + hit.start + '-' + hit.end + '\n' +
      hit.protein.match(/.{1,60}/g).join('\n') + '\n').join('')};
  }
  return {type: 'text/plain;charset=utf-8', text: '##gff-version 3\n' + hits.map((hit, index) =>
    [gffId(hit.id), 'dna-analyzer', 'open_reading_frame', hit.start, hit.end, '.', hit.frame[0], '.',
      'ID=orf' + (index + 1)].join('\t') + '\n').join('')};
}
