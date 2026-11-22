import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {test} from 'node:test';
import {orfs, parseFasta} from '../core.mjs';
import {copyRegion, mapOrf, reverseRegion, translateRegion} from '../viewer.mjs';

test('copy and reverse complement inclusive regions on both strands', () => {
  const sequence = 'CATGCCT';
  const region = {start: 2, end: 6};
  assert.equal(copyRegion(sequence, region), 'ATGCC');
  assert.equal(copyRegion(sequence, {...region, strand: '-'}), 'GGCAT');
  assert.equal(reverseRegion(sequence, region), 'GGCAT');
  assert.equal(reverseRegion(sequence, {...region, strand: '-'}), 'ATGCC');
  assert.equal(copyRegion(sequence), sequence);
  assert.equal(reverseRegion(sequence), 'AGGCATG');
});

test('translation starts at the first base of the oriented selection', () => {
  const sequence = 'CATGTGAT';
  const region = {start: 2, end: 7};
  assert.deepEqual(translateRegion(sequence, region), {frame: '+1', protein: 'M*'});
  assert.deepEqual(translateRegion(sequence, {...region, strand: '-'}), {frame: '+1', protein: 'SH'});
  assert.deepEqual(translateRegion(sequence, region, 2), {frame: '+1', protein: 'MW'});
});

test('insulin ORF arrows use inclusive reference coordinates and strand direction', () => {
  const {records} = parseFasta(readFileSync(new URL('../examples.fasta', import.meta.url), 'utf8'));
  const {sequence} = records.find(record => record.id === 'NM_000207.3');
  for (const [frame, start, end, length, from, to, direction] of [
    ['+3', 60, 392, 110, 60, 392, 'right'],
    ['-3', 206, 355, 49, 355, 206, 'left'],
  ]) {
    const hit = orfs(sequence).find(hit => hit.frame === frame && hit.start === start && hit.end === end);
    assert.equal(hit.protein.length, length);
    assert.deepEqual(mapOrf(hit), {start, end, strand: frame[0], from, to, direction});
  }
});
