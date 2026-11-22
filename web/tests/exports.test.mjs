import assert from 'node:assert/strict';
import {test} from 'node:test';
import {stats} from '../core.mjs';
import {exportSample} from '../exports.mjs';

test('exports use current ORF options and escape sample identifiers', () => {
  const record = {id: '=HYPERLINK("<test>")', sequence: 'GTGAAATAA', description: 'sample'};
  const sample = {name: '<script>', records: [record], stats: [stats(record)],
    table: 11, orfOptions: {alternative_starts: true}};
  const exportAs = format => exportSample(format, sample).text;
  assert.match(exportAs('fasta'), /\nMK\n$/);
  assert.match(exportAs('gff3'), /^##gff-version 3\n%3DHYPERLINK%28%22%3Ctest%3E%22%29\t.*\t1\t9\t/);
  assert.match(exportAs('csv'), /"'=HYPERLINK\(""<test>""\)"/);
  assert.match(exportAs('html'), /&lt;script&gt;/);
  assert.doesNotMatch(exportAs('html'), /<script>/);
});
