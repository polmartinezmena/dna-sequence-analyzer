import {codes, dna, orfs, parseFasta, stats, transcribe, translate, windows, oligoTm, OLIGO_MAX_LENGTH} from './core.mjs';
import {exportSample} from './exports.mjs';
import {blastUrl} from './blast.mjs';
import {copyRegion, mapOrf, reverseRegion, translateRegion} from './viewer.mjs';

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

const state = {name: '', records: [], stats: [], comments: 0, active: 0, page: 'viewer',
  region: null, selectedOrf: null, table: 1, orfOptions: {min_aa: 1}};
const pages = {
  viewer: ['Secuencia', 'M3 3h18v18H3zM12 3v18M3 12h18'],
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
    if (matchMedia('(max-width: 760px)').matches) {
      setMenu(false);
      $('menu-toggle').focus();
    }
    renderView();
  });
  $('navigation').append(button);
}

function table(headers, rows) {
  const head = h('thead', h('tr', headers.map(label => h('th', label, {scope: 'col'}))));
  const body = h('tbody', rows.map(values => h('tr', values.map(value => h('td', value)))));
  return h('div', h('table', [head, body]), {class: 'table-scroll'});
}

function syncRegion() {
  const length = state.stats[state.active].length;
  const {start, end} = state.region ?? {};
  const form = $('region');
  form.elements.start.value = start ?? '';
  form.elements.end.value = end ?? '';
  for (const field of [form.elements.start, form.elements.end, form.querySelector('button')]) {
    field.disabled = !length;
  }
  form.elements.start.max = length;
  form.elements.end.max = length;
  const clear = $('region-clear');
  clear.hidden = !state.region;
  $('region-status').textContent = state.region?.strand === '-' ? 'Hebra inversa' : '';
}

function applyRegion(start, end, strand = '+', selectedOrf = null) {
  state.region = {start, end, strand};
  state.selectedOrf = selectedOrf;
  state.page = 'viewer';
  syncRegion();
  renderView();
}

function pickRecord(index) {
  state.active = index;
  $('record-select').value = index;
  const length = state.stats[index].length;
  state.region = length ? {start: 1, end: Math.min(120, length), strand: '+'} : null;
  state.selectedOrf = null;
  $('sample-name').textContent = state.name;
  syncRegion();
  renderView();
}

function renderSample(record, hits, orfError) {
  const current = state.stats[state.active];
  const tooLongForTm = current.length > OLIGO_MAX_LENGTH;
  const tm = oligoTm(record.sequence);
  const metrics = [['Longitud', fmt(current.length) + ' pb'], ['GC', percent(current.gc_percent)],
    ['Tm', tooLongForTm ? 'No aplica' : tm === null ? '—' : fmt(tm) + ' °C'], ['ORFs', orfError ? '—' : fmt(hits.length)],
    ['Ambiguas', fmt(current.ambiguous_bases)]];
  $('metrics').replaceChildren(...metrics.map(([label, value]) => h('div', [h('span', label), h('strong', value)])));
  const invalid = [...new Set(record.sequence.toUpperCase().replace(/[ACGTURYSWKMBDHVN]/g, ''))];
  const checks = [
    [true, fmt(state.records.length) + ' registros FASTA leídos'],
    [Boolean(current.length), current.length ? 'Registro con secuencia' : 'Registro vacío'],
    [!invalid.length, invalid.length ? 'Fuera de IUPAC: ' + invalid.join(', ') : 'Alfabeto IUPAC válido'],
    [!current.ambiguous_bases, fmt(current.ambiguous_bases) + ' bases ambiguas'],
    [!current.gc_flag, current.gc_flag === 'low-gc' ? 'GC bajo (< 25 %)' :
      current.gc_flag === 'high-gc' ? 'GC alto (> 65 %)' : 'Sin valores extremos de GC'],
  ];
  if (state.comments) checks.push([true, fmt(state.comments) + ' comentarios de cabecera']);
  if (tooLongForTm) checks.push([true, `Tm solo para oligos de hasta ${OLIGO_MAX_LENGTH} pb`]);
  if (orfError) checks.push([false, orfError.message]);
  $('checks').replaceChildren(...checks.map(([ok, label]) => h('li', [
    h('span', ok ? '✓' : '!', {'aria-label': ok ? 'Correcto' : 'Aviso'}), label,
  ])));
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
  Object.assign(state, parsed, {name, stats: parsed.records.map(stats)});
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

function sequenceMap(record, hits) {
  const length = record.sequence.length;
  const left = 42;
  const right = 660;
  const x = position => left + (right - left) * (position - 1) / Math.max(1, length - 1);
  const trackHeight = hits.length > 2 ? 24 : 32;
  const chartTop = 64 + hits.length * trackHeight;
  const chartBottom = chartTop + 82;
  const height = chartBottom + 12;
  const map = svg('svg', {viewBox: '0 0 680 ' + height, class: 'sequence-map',
    role: 'group', 'aria-label': 'Mapa de posiciones, ORFs y GC por ventana'});
  const line = (x1, y1, x2, y2, className) => svg('line', {x1, y1, x2, y2, class: className});
  function label(value, attributes) {
    const node = svg('text', attributes);
    node.textContent = value;
    return node;
  }

  if (state.region) {
    const start = x(state.region.start);
    const width = x(state.region.end) - start;
    map.append(svg('rect', {x: start, y: 29, width: Math.max(2, width),
      height: chartBottom - 29, class: 'map-selection'}));
  }
  map.append(line(left, 29, right, 29, 'axis'));
  const interval = 10 ** Math.floor(Math.log10(length));
  const step = length / interval < 3 ? Math.max(1, Math.round(interval / 2)) : interval;
  const ticks = [1];
  for (let value = step; value < length; value += step) {
    if (value - 1 > length * .065 && length - value > length * .065) ticks.push(value);
  }
  if (length > 1) ticks.push(length);
  for (const value of ticks) {
    map.append(line(x(value), 22, x(value), 29, 'axis'),
      label(fmt(value), {x: x(value), y: 16, 'text-anchor': 'middle'}));
  }

  const principal = hits.reduce((best, hit) => !best || hit.protein.length > best.protein.length ? hit : best, null);
  hits.forEach((hit, index) => {
    const {from, to, strand, direction} = mapOrf(hit);
    const y = 54 + index * trackHeight;
    const tail = x(from);
    const tip = x(to);
    const arrowBase = tip + (direction === 'right' ? -1 : 1) * Math.min(10, Math.abs(tip - tail));
    const caption = 'ORF ' + fmt(index + 1);
    const active = state.selectedOrf
      ? ['start', 'end', 'frame'].every(key => state.selectedOrf[key] === hit[key])
      : principal === hit;
    const group = svg('g', {class: 'orf-track' + (active ? ' active' : ''),
      role: 'button', tabindex: 0, 'aria-label': caption + ', marco ' + hit.frame + ', ' + fmt(hit.start) + '–' + fmt(hit.end)});
    const hitWidth = Math.max(38, Math.abs(tail - tip));
    group.append(svg('rect', {x: (tail + tip - hitWidth) / 2, y: y - trackHeight + 5,
      width: hitWidth, height: trackHeight, fill: 'transparent'}),
      label(hit.frame, {x: 4, y: y + 4}), line(tail, y, tip, y, 'orf-line'),
      svg('polygon', {points: tip + ',' + y + ' ' + arrowBase + ',' + (y - 5) + ' ' + arrowBase + ',' + (y + 5)}),
      label(caption, {x: (tail + tip) / 2, y: y - 6, 'text-anchor': 'middle'}));
    const choose = () => {
      applyRegion(hit.start, hit.end, strand, hit);
      document.querySelectorAll('.orf-track')[index].focus();
    };
    group.addEventListener('click', event => { event.stopPropagation(); choose(); });
    group.addEventListener('keydown', event => {
      if (event.key !== 'Enter' && event.key !== ' ') return;
      event.preventDefault();
      choose();
    });
    map.append(group);
  });

  const width = Math.min(20, length);
  const rows = windows([record], width, Math.max(1, Math.ceil(length / 650)));
  const values = rows.flatMap(row => row.gc_percent === null ? [] : [row.gc_percent]);
  const minimum = Math.max(0, Math.floor(Math.min(40, ...values) / 10) * 10);
  const maximum = Math.min(100, Math.ceil(Math.max(70, ...values) / 10) * 10);
  const y = value => chartBottom - (value - minimum) * 82 / (maximum - minimum);
  for (let value = maximum; value >= minimum; value -= maximum - minimum > 50 ? 20 : 10) {
    map.append(line(left, y(value), right, y(value), 'grid'), label(fmt(value) + '%', {x: 0, y: y(value) + 4}));
  }
  map.append(line(left, chartTop - 3, left, chartBottom + 3, 'axis'));
  let connected = false;
  let path = '';
  for (const row of rows) {
    if (row.gc_percent === null) { connected = false; continue; }
    const point = x((row.start + row.end) / 2) + ',' + y(row.gc_percent);
    path += (connected ? 'L' : 'M') + point;
    if (!connected) path += 'l.01,0';
    connected = true;
  }
  map.append(svg('path', {d: path, class: 'gc-line'}));
  map.addEventListener('click', event => {
    const at = new DOMPoint(event.clientX, event.clientY).matrixTransform(map.getScreenCTM().inverse()).x;
    const position = Math.max(1, Math.min(length, 1 + Math.round((at - left) * (length - 1) / (right - left))));
    applyRegion(position, position);
  });
  return map;
}

function viewer(record, hits, firstBase) {
  const card = h('div', undefined, {class: 'card viewer'});
  const length = record.sequence.length;
  if (!length) return h('div', 'El registro está vacío.', {class: 'card padded'});
  card.append(h('div', [h('h2', 'Mapa'), sequenceMap(record, hits)],
    {class: 'map-panel', tabindex: 0, 'aria-label': 'Mapa de la secuencia'}));

  const output = h('div', undefined, {class: 'action-output', role: 'status'});
  const toolbar = h('div', undefined, {class: 'toolbar', role: 'group', 'aria-label': 'Acciones de secuencia'});
  const result = (label, value) => {
    const field = h('textarea', undefined, {rows: 2, readonly: '', class: 'mono'});
    field.value = value;
    output.replaceChildren(h('label', [label, field]));
  };
  const actions = [
    ['Copiar', async () => {
      const value = copyRegion(record.sequence, state.region);
      try {
        await navigator.clipboard.writeText(value);
        output.textContent = 'Secuencia copiada.';
      } catch {
        result('Secuencia para copiar', value);
        output.querySelector('textarea').select();
      }
    }],
    ['Complemento inverso', () => result('Complemento inverso', reverseRegion(record.sequence, state.region))],
    ['Traducir', () => {
      const {protein} = translateRegion(record.sequence, state.region, state.table);
      result('Traducción desde la primera base de la hebra activa · tabla ' + fmt(state.table), protein);
    }],
  ];
  for (const [label, action] of actions) {
    const button = h('button', label, {type: 'button'});
    button.addEventListener('click', async () => {
      try { await action(); }
      catch (error) { output.textContent = error.message; }
    });
    toolbar.append(button);
  }

  const scroll = h('div', undefined, {class: 'sequence', tabindex: 0, 'aria-label': 'Visor de bases'});
  const columns = matchMedia('(max-width: 540px)').matches ? 24
    : matchMedia('(max-width: 1120px)').matches ? 36 : 60;
  // Bound DOM size while keeping the selected interval reachable through the map and form.
  const visibleBases = columns * 24;
  const first = firstBase ?? Math.floor(((state.region?.start ?? 1) - 1) / visibleBases) * visibleBases;
  const last = Math.min(length, first + visibleBases);
  let selectedRow;
  for (let pos = first; pos < last; pos += columns) {
    const bases = h('span', undefined, {class: 'codons'});
    for (let i = pos; i < Math.min(pos + columns, last); i += 3) {
      const codon = h('span', undefined, {class: 'codon'});
      for (let j = i; j < Math.min(i + 3, last); j++) {
        const selected = state.region && j + 1 >= state.region.start && j + 1 <= state.region.end;
        codon.append(h('span', record.sequence[j].toUpperCase(), {class: selected ? 'selected-base' : ''}));
      }
      bases.append(codon);
    }
    const row = h('div', [h('span', fmt(pos + 1), {class: 'position'}), bases], {class: 'sequence-row'});
    if (state.region && pos < state.region.start && state.region.start <= pos + columns) selectedRow = row;
    scroll.append(row);
  }
  const viewerPanel = h('div', [toolbar, output, scroll], {class: 'viewer-panel'});
  if (first || last < length) {
    const range = h('div', h('span', 'Bases ' + fmt(first + 1) + '–' + fmt(last) + ' de ' + fmt(length)),
      {class: 'sequence-range'});
    for (const [label, offset, unavailable] of [
      ['Anteriores', Math.max(0, first - visibleBases), first === 0],
      ['Siguientes', last, last === length],
    ]) {
      const button = h('button', label, {type: 'button'});
      button.disabled = unavailable;
      button.addEventListener('click', () => {
        const next = viewer(record, hits, offset);
        card.replaceWith(next);
        next.querySelector('.sequence').focus();
      });
      range.append(button);
    }
    viewerPanel.append(range);
  }
  card.append(viewerPanel);
  if (selectedRow) requestAnimationFrame(() => {
    if (selectedRow.offsetTop + selectedRow.offsetHeight > scroll.clientHeight) {
      scroll.scrollTop = selectedRow.offsetTop;
    }
  });
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
  let hits = [];
  let orfError;
  try { hits = getOrfs(record); }
  catch (error) { orfError = error; }
  renderSample(record, hits, orfError);
  if (state.page === 'viewer') return content.replaceChildren(viewer(record, hits));
  const card = h('div', h('h2', pages[state.page][0]), {class: 'card padded section-content'});
  content.replaceChildren(card);
  try {

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
  if (end < start) {
    $('region-status').textContent = 'El fin debe ser igual o posterior al inicio.';
    return;
  }
  applyRegion(start, end);
});
$('region-clear').addEventListener('click', () => {
  state.region = null;
  state.selectedOrf = null;
  syncRegion();
  renderView();
  $('region').elements.start.focus();
});
for (const rule of ['(max-width: 540px)', '(max-width: 1120px)']) {
  matchMedia(rule).addEventListener('change', () => {
    if (state.page === 'viewer') renderView();
  });
}

function setMenu(open) {
  $('menu-toggle').setAttribute('aria-expanded', open);
  $('sidebar-controls').classList.toggle('open', open);
}
$('menu-toggle').addEventListener('click', () => {
  setMenu($('menu-toggle').getAttribute('aria-expanded') !== 'true');
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
