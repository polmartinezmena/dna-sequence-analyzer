import {reverseComplement, translate} from './core.mjs';

export function copyRegion(sequence, region) {
  if (!region) return sequence;
  const bases = sequence.slice(region.start - 1, region.end);
  return region.strand === '-' ? reverseComplement(bases) : bases;
}

export function reverseRegion(sequence, region) {
  const bases = region ? sequence.slice(region.start - 1, region.end) : sequence;
  return region?.strand === '-' ? bases : reverseComplement(bases);
}

export function translateRegion(sequence, region, table = 1) {
  return translate(copyRegion(sequence, region), table, '+1')[0];
}

export function mapOrf({start, end, frame}) {
  const strand = frame[0];
  return strand === '+'
    ? {start, end, strand, from: start, to: end, direction: 'right'}
    : {start, end, strand, from: end, to: start, direction: 'left'};
}
