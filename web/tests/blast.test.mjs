import assert from 'node:assert/strict';
import {test} from 'node:test';
import {blastUrl} from '../blast.mjs';

test('NCBI opens a prefilled form without submitting a job', () => {
  const url = new URL(blastUrl('auGC'));
  assert.equal(url.origin, 'https://blast.ncbi.nlm.nih.gov');
  assert.equal(url.searchParams.get('PAGE_TYPE'), 'BlastSearch');
  assert.equal(url.searchParams.get('PROGRAM'), 'blastn');
  assert.equal(url.searchParams.get('QUERY'), 'ATGC');
  assert.equal(url.searchParams.get('MEGABLAST'), 'on');
  assert.equal(url.searchParams.has('CMD'), false);
  const orf = new URL(blastUrl('MKL', 'blastp'));
  assert.equal(orf.searchParams.get('PROGRAM'), 'blastp');
  assert.equal(orf.searchParams.get('QUERY'), 'MKL');
  assert.equal(orf.searchParams.has('MEGABLAST'), false);
});

test('long queries use manual paste and invalid DNA fails', () => {
  assert.equal(new URL(blastUrl('A'.repeat(2000))).searchParams.has('QUERY'), false);
  assert.throws(() => blastUrl('ACG-'));
});
