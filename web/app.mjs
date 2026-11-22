import {codes, dna, orfs, parseFasta, reverseComplement, stats, transcribe, translate, windows, melt} from './core.mjs';
import {exportSample} from './exports.mjs';
import {blastUrl} from './blast.mjs';

const $ = id => document.getElementById(id);
const h = (tag, content, attrs = {}) => {
  const element = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
  if (content !== undefined) element.append(...[].concat(content));
  return element;
};
const number = new Intl.NumberFormat('es-ES', {maximumFractionDigits: 2});
const fmt = value => value == null ? '—' : number.format(value);
const percent = value => value == null ? '—' : fmt(value) + ' %';
const regionLimit = 500;
const state = {name: '', records: [], stats: [], comments: 0, active: 0, page: 'viewer',
  tab: 0, start: 1, end: 1, table: 1, orfOptions: {min_aa: 1}};
const pages = {
  summary: ['Resumen', 'M3 3h8v8H3zM14 3h7v8h-7zM3 14h8v7H3zM14 14h7v7h-7z'],
  viewer: ['Visor', 'M3 5h18v14H3zM7 5v14m5-14v14m5-14v14'],
  translation: ['Traducción', 'M3 7h8m-4-4v4m3 0c0 4-3 7-6 9m1-7c1 2 3 4 6 5m1 7 4-10 4 10m-7-3h6'],
  orfs: ['ORFs', 'M3 6h18M3 12h18M3 18h18M7 4v4m9 2v4m-6 2v4'],
  blast: ['BLAST', 'M16 16l5 5M18 10a8 8 0 1 1-16 0 8 8 0 0 1 16 0Z'],
};
const svgNS = 'http://www.w3.org/2000/svg';
function svg(tag, attributes) {
  const element = document.createElementNS(svgNS, tag);
  for (const [key, value] of Object.entries(attributes)) element.setAttribute(key, value);
  return element;
}

for (const [page, [label, path]] of Object.entries(pages)) {
  const icon = svg('svg', {viewBox: '0 0 24 24', 'aria-hidden': 'true'});
  icon.append(svg('path', {d: path}));
  const button = h('button', [icon, label], {type: 'button', id: 'nav-' + page});
  button.addEventListener('click', () => {
    state.page = page;
    renderView();
  });
  $('navigation').append(button);
}

function table(headers, rows) {
  const head = h('thead', h('tr', headers.map(label => h('th', label, {scope: 'col'}))));
  const body = h('tbody', rows.map(values => h('tr', values.map(value => h('td', value)))));
  return h('div', h('table', [head, body]), {class: 'table-scroll'});
}

function resetRegion() {
  const length = state.stats[state.active].length;
  state.start = 1;
  state.end = Math.min(120, length);
  const form = $('region');
  form.elements.start.value = length ? state.start : '';
  form.elements.end.value = length ? state.end : '';
  for (const field of [form.elements.start, form.elements.end, form.querySelector('button')]) {
    field.disabled = !length;
  }
  form.elements.start.max = length;
  form.elements.end.max = length;
  $('region-status').textContent = '';
}

function pickRecord(index) {
  state.active = index;
  $('record-select').value = index;
  resetRegion();
  renderSample();
  renderView();
}

function renderSample() {
  const total = key => state.stats.reduce((sum, item) => sum + item[key], 0);
  const countable = total('unambiguous_bases');
  const metrics = [
    ['Registros', fmt(state.records.length)], ['Bases totales', fmt(total('length'))],
    ['GC medio', countable ? percent(100 * total('gc_bases') / countable) : '—'],
  ];
  $('filename').textContent = $('sample-name').textContent = state.name;
  $('metrics').replaceChildren(...metrics.map(([label, value]) => h('div', [h('span', label), h('strong', value)])));
  const current = state.stats[state.active];
  const invalid = [...new Set(state.records[state.active].sequence.toUpperCase().replace(/[ACGTURYSWKMBDHVN]/g, ''))];
  const checks = [
    [true, fmt(state.records.length) + ' registros FASTA leídos'],
    [Boolean(current.length), current.length ? 'Registro con secuencia' : 'Registro vacío'],
    [!invalid.length, invalid.length ? 'Fuera de IUPAC: ' + invalid.join(', ') : 'Alfabeto IUPAC válido'],
    [!current.ambiguous_bases, fmt(current.ambiguous_bases) + ' bases ambiguas'],
    [!current.gc_flag, current.gc_flag === 'low-gc' ? 'GC bajo (< 25 %)' :
      current.gc_flag === 'high-gc' ? 'GC alto (> 65 %)' : 'Sin valores extremos de GC'],
  ];
  if (state.comments) checks.push([true, fmt(state.comments) + ' comentarios de cabecera']);
  $('checks').replaceChildren(...checks.map(([ok, label]) => h('li', [
    h('span', ok ? '✓' : '!', {class: ok ? '' : 'warn', 'aria-label': ok ? 'Correcto' : 'Aviso'}), label,
  ])));
  const details = [['Identificador', current.id || '(sin ID)'], ['Longitud', fmt(current.length) + ' pb'],
    ['GC', percent(current.gc_percent)], ['Ambiguas', fmt(current.ambiguous_bases)]];
  $('record-details').replaceChildren(...details.map(([key, value]) => h('div', [h('dt', key), h('dd', value)])));

  const rows = state.stats.map((item, i) => {
    const button = h('button', item.id || '(sin ID)', {type: 'button', 'aria-pressed': i === state.active});
    button.addEventListener('click', () => {
      pickRecord(i);
      $('records').querySelectorAll('button')[i].focus();
    });
    return [button, fmt(item.length), percent(item.gc_percent), fmt(item.ambiguous_bases)];
  });
  const grid = table(['Identificador', 'Longitud (pb)', 'GC', 'Bases ambiguas'], rows);
  grid.querySelectorAll('tbody tr')[state.active].setAttribute('data-selected', '');
  $('records').replaceChildren(grid);
}

function load(text, name, raw = false) {
  let parsed;
  const stripped = text.trimStart();
  if (raw && !/^(>|[;#])/.test(stripped)) {
    const sequence = stripped.replace(/\s/g, '');
    if (!/[ACGTU]/i.test(sequence)) throw new Error('Pega un FASTA o una secuencia de bases.');
    parsed = {comments: 0, records: [{id: 'secuencia', description: 'secuencia', sequence}]};
  } else {
    parsed = parseFasta(text);
  }
  if (!parsed.records.length) throw new Error('El archivo no contiene registros FASTA.');
  Object.assign(state, parsed, {name, stats: parsed.records.map(stats), tab: 0});
  const totals = new Map();
  for (const {id} of parsed.records) totals.set(id, (totals.get(id) ?? 0) + 1);
  const seen = new Map();
  $('record-select').replaceChildren(...parsed.records.map(({id}, index) => {
    const count = (seen.get(id) ?? 0) + 1;
    seen.set(id, count);
    return h('option', (id || '(sin ID)') + (totals.get(id) > 1 ? ' (' + count + ')' : ''), {value: index});
  }));
  $('status').textContent = '';
  $('export-open').disabled = false;
  pickRecord(0);
}

function chart(rows, key, label, range) {
  const values = rows.map(row => row[key]).filter(value => value !== null);
  if (!values.length) return h('p', 'Sin datos para ' + label + '.');
  const [min, max] = range ?? [Math.min(0, ...values), Math.max(1, ...values)];
  const height = 113;
  const element = svg('svg', {viewBox: '0 0 600 132', class: 'chart', role: 'img', 'aria-label': label});
  const first = rows[0].start;
  const span = Math.max(1, rows.at(-1).start - first);
  for (const [value, y] of [[max, 13], [(min + max) / 2, 57], [min, 101]]) {
    element.append(svg('line', {x1: 44, x2: 574, y1: y, y2: y, class: 'grid'}));
    const tick = svg('text', {x: 1, y: y + 4});
    tick.textContent = fmt(value);
    element.append(tick);
  }
  let path = '';
  let started = false;
  for (const row of rows) {
    const value = row[key];
    if (value === null) { started = false; continue; }
    const x = 44 + (row.start - first) * 530 / span;
    const y = 101 - (value - min) * 88 / (max - min);
    path += (started ? 'L' : 'M') + x + ',' + y;
    if (!started) element.append(svg('circle', {cx: x, cy: y, r: 2, fill: 'currentColor'}));
    started = true;
  }
  element.append(svg('path', {d: path, class: 'line'}));
  for (const [position, x, anchor] of [[first, 44, 'start'], [rows.at(-1).start, 574, 'end']]) {
    const tick = svg('text', {x, y: height + 12, 'text-anchor': anchor});
    tick.textContent = fmt(position);
    element.append(tick);
  }
  return element;
}

function viewer(record) {
  const card = h('div', undefined, {class: 'viewer'});
  const tablist = h('div', undefined, {class: 'tabs', role: 'tablist', 'aria-label': 'Vistas de secuencia'});
  const tabs = ['Secuencia', 'Vista por tramos', 'Anotaciones'].map((label, index) => {
    const button = h('button', label, {id: 'tab-' + index, role: 'tab', type: 'button',
      'aria-selected': index === state.tab, 'aria-controls': 'viewer-panel', tabindex: index === state.tab ? 0 : -1});
    button.addEventListener('click', () => {
      state.tab = index;
      renderView();
      $('tab-' + index).focus();
    });
    button.addEventListener('keydown', event => {
      const target = {ArrowLeft: (index + 2) % 3, ArrowRight: (index + 1) % 3, Home: 0, End: 2}[event.key];
      if (target === undefined) return;
      event.preventDefault();
      tabs[target].click();
    });
    return button;
  });
  tablist.append(...tabs);
  const panel = h('div', undefined, {id: 'viewer-panel', role: 'tabpanel', tabindex: '0',
    'aria-labelledby': 'tab-' + state.tab, class: 'section-content'});
  card.append(h('h2', 'Visor de secuencia'), tablist, panel);
  if (!record.sequence.length) {
    panel.append(h('p', 'El registro está vacío.'));
    return card;
  }
  const sequence = record.sequence.slice(state.start - 1, state.end).toUpperCase();
  panel.append(h('p', record.id + ' · ' + fmt(state.start) + '–' + fmt(state.end) + ' pb'));
  try {
    if (state.tab === 2) {
      const hits = getOrfs(record).filter(hit => hit.start <= state.end && hit.end >= state.start);
      panel.append(orfTable(hits));
    } else if (state.tab === 1) {
      panel.append(h('h3', 'Hebra directa 5′ → 3′'), h('pre', sequence, {class: 'mono'}),
        h('h3', 'Complementaria inversa 5′ → 3′'), h('pre', reverseComplement(sequence), {class: 'mono'}));
    } else {
      const lines = h('div', undefined, {class: 'sequence', 'aria-label': 'Bases de la región'});
      const columns = matchMedia('(max-width: 540px)').matches ? 12 : 24;
      for (let pos = 0; pos < sequence.length; pos += columns) {
        const start = state.start + pos;
        const ruler = h('div', undefined, {class: 'lane ruler', style: '--columns:' + columns, 'aria-hidden': 'true'});
        const lane = h('div', undefined, {class: 'lane', style: '--columns:' + columns});
        for (let col = 0; col < columns && pos + col < sequence.length; col++) {
          const base = sequence[pos + col];
          ruler.append(h('span', col === 0 || (start + col) % 10 === 0 ? String(start + col) : ''));
          lane.append(h('span', base, {class: 'tile', 'data-base': base, title: 'Base ' + fmt(start + col)}));
        }
        lines.append(ruler, lane);
      }
      const width = Math.min(20, sequence.length);
      const rows = windows([{id: record.id, sequence}], width, 1);
      for (const row of rows) row.start += state.start - 1;
      const label = 'GC por posición (%) · ventana de ' + fmt(width) + ' pb';
      panel.append(lines, h('h3', label), chart(rows, 'gc_percent', label, [0, 100]));
    }
  } catch (error) {
    panel.append(h('p', error.message, {role: 'status'}));
  }
  return card;
}

function getOrfs(record) {
  return orfs(record.sequence, {...state.orfOptions, table: state.table});
}

function orfTable(hits) {
  if (!hits.length) return h('p', 'Sin ORFs con estos parámetros.');
  return table(['Marco', 'Inicio', 'Fin', 'Longitud (aa)', 'Codón inicial', 'Parada'], hits.map(hit => [
    hit.frame, fmt(hit.start), fmt(hit.end), fmt(hit.protein.length), hit.start_codon, hit.has_stop ? 'Sí' : 'Parcial',
  ]));
}

function tableSelect() {
  const select = h('select', Object.entries(codes).map(([id, code]) => h('option', id + ' · ' + code.name, {value: id})));
  select.value = state.table;
  select.addEventListener('change', () => { state.table = Number(select.value); renderView(); $('content').querySelector('select').focus(); });
  return h('label', ['Tabla genética NCBI', select]);
}

function summary(record, card) {
  const length = record.sequence.length;
  if (!length) return card.append(h('p', 'El registro está vacío.'));
  const width = Math.min(60, length);
  const step = Math.max(1, Math.ceil(length / 600));
  const rows = windows([record], width, step);
  card.append(h('p', 'Ventana de ' + fmt(width) + ' pb · Paso de ' + fmt(step) + ' pb'));
  for (const [key, label, bounds] of [['gc_percent', 'GC por posición (%)', [0, 100]],
    ['cpg_oe', 'CpG observado / esperado'], ['gc_skew', 'Sesgo GC', [-1, 1]],
    ['cumulative_skew', 'Sesgo GC acumulado']]) {
    card.append(h('h3', label), chart(rows, key, label, bounds));
  }
}

function blast(record, card) {
  card.append(h('p', 'NCBI abre el formulario con tu consulta. Revisa los parámetros y pulsa BLAST allí. El navegador envía la secuencia a NCBI al abrir el enlace.'));
  const sources = new Map([['record', {sequence: record.sequence, program: 'blastn'}]]);
  const options = [h('option', 'Registro completo · blastn', {value: 'record'})];
  try {
    getOrfs(record).forEach((hit, index) => {
      const key = String(index);
      sources.set(key, {sequence: hit.protein, program: 'blastp'});
      options.push(h('option', 'ORF ' + hit.frame + ' · ' + fmt(hit.start) + '–' + fmt(hit.end) +
        ' · ' + fmt(hit.protein.length) + ' aa', {value: key}));
    });
  } catch (error) {
    card.append(h('p', 'ORFs no disponibles: ' + error.message));
  }
  const select = h('select', options);
  const action = h('div', undefined, {class: 'section-content'});
  function update() {
    const {sequence, program} = sources.get(select.value);
    action.replaceChildren();
    try {
      const url = new URL(blastUrl(sequence, program));
      if (!url.searchParams.has('QUERY')) {
        const query = h('textarea', undefined, {rows: 4, readonly: '', class: 'mono'});
        query.value = program === 'blastn' ? dna(sequence) : sequence;
        query.addEventListener('focus', () => query.select());
        action.append(h('p', 'La consulta supera el límite del enlace. Cópiala y pégala en NCBI.'),
          h('label', ['Secuencia para copiar', query]));
      }
      action.append(h('a', 'Abrir NCBI BLAST ↗', {href: url.href, target: '_blank',
        rel: 'noopener noreferrer', class: 'button primary blast-link'}));
    } catch (error) {
      action.append(h('p', error.message, {role: 'status'}));
    }
  }
  select.addEventListener('change', update);
  card.append(h('label', ['Consulta', select]), action);
  update();
}


function renderView() {
  for (const page of Object.keys(pages)) {
    const button = $('nav-' + page);
    if (page === state.page) button.setAttribute('aria-current', 'page');
    else button.removeAttribute('aria-current');
  }
  const record = state.records[state.active];
  const content = $('content');
  if (!record) return;
  if (state.page === 'viewer') return content.replaceChildren(viewer(record));
  const card = h('div', h('h2', pages[state.page][0]), {class: 'card section-content'});
  content.replaceChildren(card);
  try {
    if (state.page === 'summary') summary(record, card);
    if (state.page === 'blast') blast(record, card);
    if (state.page === 'translation') {
      card.append(tableSelect(), h('h3', 'ARN transcrito'), h('pre', transcribe(record.sequence), {class: 'mono'}));
      for (const frame of translate(record.sequence, state.table)) {
        card.append(h('h3', 'Marco ' + frame.frame + ' · ' + fmt(frame.protein.length) + ' aa'),
          h('pre', frame.protein, {class: 'mono'}));
      }
    }
    if (state.page === 'orfs') {
      const form = h('form', [
        h('label', ['Mínimo de aminoácidos', h('input', undefined, {name: 'min_aa', type: 'number',
          min: 1, value: state.orfOptions.min_aa, required: ''})]),
        ...Object.entries({alternative_starts: 'Inicios alternativos', include_nested: 'Incluir anidados',
          allow_no_stop: 'Incluir parciales'}).map(([name, label]) => {
          const input = h('input', undefined, {type: 'checkbox', name});
          input.checked = Boolean(state.orfOptions[name]);
          return h('label', [input, ' ' + label]);
        }),
        h('button', 'Buscar ORFs', {type: 'submit'}),
      ], {class: 'section-content'});
      form.addEventListener('submit', event => {
        event.preventDefault();
        state.orfOptions = {
          min_aa: form.elements.min_aa.valueAsNumber, alternative_starts: form.elements.alternative_starts.checked,
          include_nested: form.elements.include_nested.checked, allow_no_stop: form.elements.allow_no_stop.checked,
        };
        renderView();
        $('content').querySelector('form button').focus();
      });
      card.append(tableSelect(), form, orfTable(getOrfs(record)));
    }

  } catch (error) {
    card.append(h('p', error.message, {role: 'status'}));
  }
}

$('record-select').addEventListener('change', event => pickRecord(Number(event.target.value)));
$('region').addEventListener('submit', event => {
  event.preventDefault();
  const start = event.target.elements.start.valueAsNumber;
  const end = event.target.elements.end.valueAsNumber;
  if (end < start || end - start + 1 > regionLimit) {
    $('region-status').textContent = 'Elige una región de 1 a 500 bases con el fin igual o posterior al inicio.';
    return;
  }
  Object.assign(state, {start, end, page: 'viewer'});
  $('region-status').textContent = 'Región aplicada.';
  renderView();
});
matchMedia('(max-width: 540px)').addEventListener('change', () => {
  if (state.page === 'viewer') renderView();
});

const theme = $('theme');
function setTheme(dark) {
  document.documentElement.dataset.theme = dark ? 'dark' : 'light';
  theme.setAttribute('aria-pressed', dark);
  theme.textContent = dark ? 'Tema claro' : 'Tema oscuro';
  document.querySelector('meta[name="theme-color"]').content = dark ? '#171411' : '#F6F4EF';
}
let savedTheme;
try { savedTheme = localStorage.getItem('dna-theme'); } catch {}
setTheme(savedTheme ? savedTheme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches);
theme.addEventListener('click', () => {
  const dark = document.documentElement.dataset.theme !== 'dark';
  setTheme(dark);
  try { localStorage.setItem('dna-theme', dark ? 'dark' : 'light'); } catch {}
});

const importDialog = $('import-dialog');
$('import-open').addEventListener('click', () => {
  $('import-status').textContent = '';
  importDialog.showModal();
});
function receive(text, name, raw = false) {
  try {
    load(text, name, raw);
    if (importDialog.open) importDialog.close();
  } catch (error) {
    (importDialog.open ? $('import-status') : $('status')).textContent = error.message;
  }
}
function readFile(file) {
  if (!file) return;
  const reader = new FileReader();
  reader.addEventListener('load', () => receive(reader.result, file.name));
  reader.addEventListener('error', () => {
    (importDialog.open ? $('import-status') : $('status')).textContent = 'No se pudo leer el archivo.';
  });
  reader.readAsText(file);
}
$('file').addEventListener('change', event => {
  readFile(event.target.files[0]);
  event.target.value = '';
});
document.addEventListener('dragover', event => {
  if (event.dataTransfer.types.includes('Files')) event.preventDefault();
});
document.addEventListener('drop', event => {
  if (!event.dataTransfer.files.length) return;
  event.preventDefault();
  readFile(event.dataTransfer.files[0]);
});
$('paste-form').addEventListener('submit', event => {
  event.preventDefault();
  receive($('paste').value, 'pegado.fasta', true);
});
async function example() {
  try {
    const response = await fetch('./examples.fasta');
    if (!response.ok) throw new Error('No se pudo cargar el ejemplo.');
    const text = await response.text();
    receive(text, 'examples.fasta');
  } catch {
    (importDialog.open ? $('import-status') : $('status')).textContent = 'No se pudo cargar el ejemplo.';
  }
}
$('example').addEventListener('click', example);
const exportDialog = $('export-dialog');
$('export-open').disabled = true;
$('export-open').addEventListener('click', () => {
  $('export-status').textContent = '';
  exportDialog.showModal();
});
exportDialog.addEventListener('click', event => {
  const format = event.target.dataset.format;
  if (!format) return;
  try {
    const payload = exportSample(format, state);
    const url = URL.createObjectURL(new Blob([payload.text], {type: payload.type}));
    const link = h('a', undefined, {href: url, download: 'dna-analysis.' + format});
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    $('export-status').textContent = 'Descarga preparada.';
  } catch (error) { $('export-status').textContent = error.message; }
});

renderView();
example();
if ('serviceWorker' in navigator) {
  navigator.serviceWorker.register('./sw.js').catch(error => console.warn('Service worker:', error.message));
}
