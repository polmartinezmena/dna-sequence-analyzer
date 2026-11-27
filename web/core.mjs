import {codes, complement, iupac} from './reference.mjs';
import {melt} from './melting.mjs';
export {codes, melt};

// Same limit as OLIGO_GUIDANCE_LENGTH in tm.py: the nearest-neighbor model is for short oligos.
export const OLIGO_MAX_LENGTH = 60;

// Tm of a short oligo, or null when the sequence is too long for the model or cannot be scored.
export function oligoTm(sequence) {
  return sequence.length && sequence.length <= OLIGO_MAX_LENGTH ? melt(sequence) : null;
}

export function parseFasta(text) {
  const records = [];
  let comments = 0;
  for (const line of text.replace(/^\uFEFF/, '').split(/\r\n|\n|\r/)) {
    if (line.startsWith('>')) {
      const description = line.slice(1).trimEnd();
      records.push({id: description.trim().split(/\s+/)[0], description, chunks: []});
    } else if (records.length) {
      records.at(-1).chunks.push(line.trimEnd().replace(/[ \t]/g, ''));
    } else if (/^[;#]/.test(line.trimStart())) {
      comments++;
    } else if (line.trim()) {
      throw new Error("Cada registro FASTA debe empezar por '>'.");
    }
  }
  return {comments, records: records.map(({chunks, ...record}) => ({...record, sequence: chunks.join('')}))};
}

export function stats({id, description, sequence}) {
  let gc_bases = 0;
  let at = 0;
  let ambiguous_bases = 0;
  const upper = sequence.toUpperCase();
  for (const base of upper) {
    if ('GCS'.includes(base)) gc_bases++;
    else if ('ATUW'.includes(base)) at++;
    else if ('RYKMBDHVNX'.includes(base)) ambiguous_bases++;
  }
  const unambiguous_bases = gc_bases + at;
  const gc_percent = unambiguous_bases ? 100 * gc_bases / unambiguous_bases : null;
  const hasT = upper.includes('T');
  const hasU = upper.includes('U');
  const nucleic_acid = hasT ? (hasU ? 'mixed' : 'dna') : (hasU ? 'rna' : 'unspecified');
  let gc_flag = null;
  if (gc_percent !== null && gc_percent < 25) gc_flag = 'low-gc';
  if (gc_percent > 65) gc_flag = 'high-gc';
  return {id, description, length: sequence.length, gc_bases, unambiguous_bases,
    ambiguous_bases, gc_percent, nucleic_acid, gc_flag};
}

export function windows(records, window, step) {
  for (const span of [window, step]) {
    if (!Number.isInteger(span) || span < 1) throw new Error('Usa enteros positivos para ventana y paso.');
  }
  const rows = [];
  for (const {id, sequence} of records) {
    let cumulative_skew = null;
    for (let start = 0; start + window <= sequence.length; start += step) {
      const end = start + window;
      const chunk = sequence.slice(start, end);
      const {gc_percent, ambiguous_bases} = stats({sequence: chunk});
      const clean = chunk.replace(/\s/g, '').toUpperCase();
      const c = (clean.match(/C/g) ?? []).length;
      const g = (clean.match(/G/g) ?? []).length;
      const cg = (clean.match(/CG/g) ?? []).length;
      const cpg_oe = c && g ? cg * clean.length / (c * g) : null;
      const gc_skew = c + g ? (g - c) / (g + c) : null;
      if (gc_skew !== null) cumulative_skew = (cumulative_skew ?? 0) + gc_skew;
      rows.push({id, start: start + 1, end, gc_percent, ambiguous_bases, cpg_oe, gc_skew, cumulative_skew});
    }
  }
  return rows;
}

export function dna(sequence) {
  const normalized = sequence.toUpperCase().replaceAll('U', 'T');
  const invalid = [...new Set([...normalized].filter(base => !/^[ACGTRYSWKMBDHVN]$/.test(base)))].sort();
  if (invalid.length) throw new Error('Bases desconocidas: ' + invalid.join(', '));
  return normalized;
}

export const transcribe = sequence => dna(sequence).replaceAll('T', 'U');

function reverse(dnaSequence) {
  return Array.from(dnaSequence, base => complement[base]).reverse().join('');
}

export const reverseComplement = sequence => reverse(dna(sequence));
const frames = ['+1', '+2', '+3', '-1', '-2', '-3'];
const codonCache = new WeakMap();

function geneticCode(table) {
  if (!Object.hasOwn(codes, table)) throw new Error('Tabla NCBI desconocida: ' + table);
  return codes[table];
}

function protein(sequence, offset, end, code) {
  if (!codonCache.has(code)) codonCache.set(code, new Map(Object.entries(code.amino)));
  const cache = codonCache.get(code);
  const residues = [];
  for (let pos = offset; pos + 3 <= end; pos += 3) {
    const triplet = sequence.slice(pos, pos + 3);
    if (!cache.has(triplet)) {
      const residuesForCodon = new Set();
      for (const a of iupac[triplet[0]]) {
        for (const b of iupac[triplet[1]]) {
          for (const c of iupac[triplet[2]]) residuesForCodon.add(code.amino[a + b + c] ?? '*');
        }
      }
      const symbols = [...residuesForCodon].sort().join('');
      cache.set(triplet, symbols.length === 1 ? symbols : ({DN: 'B', EQ: 'Z', IL: 'J'}[symbols] ?? 'X'));
    }
    residues.push(cache.get(triplet));
  }
  return residues.join('');
}

export function translate(sequence, table = 1, frame = 'all') {
  const code = geneticCode(table);
  if (frame !== 'all' && !frames.includes(frame)) throw new Error('Marco desconocido: ' + frame);
  const forward = dna(sequence);
  const backward = reverse(forward);
  return (frame === 'all' ? frames : [frame]).map(label => ({
    frame: label,
    protein: protein(label.startsWith('+') ? forward : backward, Number(label[1]) - 1, forward.length, code),
  }));
}

export function orfs(sequence, {table = 1, min_aa = 1, alternative_starts = false,
  include_nested = false, allow_no_stop = false} = {}) {
  if (!Number.isInteger(min_aa) || min_aa < 1) throw new Error('Indica al menos 1 aminoácido.');
  const code = geneticCode(table);
  const stops = new Set(code.stops.filter(codon => !Object.hasOwn(code.amino, codon)));
  if (!stops.size) throw new Error('La tabla no tiene codones de parada inequívocos.');
  const starts = new Set(alternative_starts ? code.starts : ['ATG']);
  const forward = dna(sequence);
  const backward = reverse(forward);
  const size = forward.length;
  const hits = [];
  for (const frame of frames) {
    const positive = frame.startsWith('+');
    const strand = positive ? forward : backward;
    const pending = [];
    const offset = Number(frame[1]) - 1;
    const limit = offset + 3 * Math.max(0, Math.floor((size - offset) / 3));
    for (let pos = offset; pos < limit + Number(allow_no_stop); pos += 3) {
      const has_stop = pos < limit;
      const codon = strand.slice(pos, pos + 3);
      if (!has_stop || stops.has(codon)) {
        for (const start of pending) {
          if (pos - start < 3 * min_aa) break;
          const end = has_stop ? pos + 3 : pos;
          hits.push({frame, start: (positive ? start : size - end) + 1, end: positive ? end : size - start,
            protein: 'M' + protein(strand, start + 3, pos, code), has_stop,
            start_codon: strand.slice(start, start + 3)});
        }
        pending.length = 0;
      } else if (starts.has(codon) && (include_nested || !pending.length)) {
        pending.push(pos);
      }
    }
  }
  return hits;
}
