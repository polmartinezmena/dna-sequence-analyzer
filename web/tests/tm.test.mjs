import assert from 'node:assert/strict';
import {test} from 'node:test';
import {melt, oligoTm, OLIGO_MAX_LENGTH} from '../core.mjs';

test('oligoTm matches melt for oligos up to the guidance length', () => {
  const primer = 'CGTTCCAAAGATGTGGGCATGAGCTTAC';
  assert.equal(oligoTm(primer), melt(primer));
  const edge = 'ACGT'.repeat(OLIGO_MAX_LENGTH / 4);
  assert.equal(edge.length, OLIGO_MAX_LENGTH);
  assert.equal(oligoTm(edge), melt(edge));
});

test('oligoTm does not report a Tm for sequences longer than an oligo', () => {
  assert.equal(oligoTm('ACGT'.repeat(OLIGO_MAX_LENGTH / 4) + 'A'), null);
  assert.equal(oligoTm('ACGT'.repeat(120)), null);
});

test('oligoTm returns null for an empty sequence', () => {
  assert.equal(oligoTm(''), null);
});
