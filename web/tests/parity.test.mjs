import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {test} from 'node:test';
import {parseFasta, stats, windows, transcribe, translate, orfs, melt} from '../core.mjs';

const fixtures = JSON.parse(readFileSync(new URL('./fixtures/parity.json', import.meta.url)));
function compare(actual, expected) {
  if (Array.isArray(expected)) {
    assert.ok(Array.isArray(actual));
    assert.equal(actual.length, expected.length);
    expected.forEach((value, index) => compare(actual[index], value));
  } else if (expected !== null && typeof expected === 'object') {
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort());
    for (const key in expected) {
      if (['gc_percent', 'tm'].includes(key) && expected[key] !== null) {
        assert.ok(Math.abs(actual[key] - expected[key]) <= 0.01 && actual[key] !== null);
      } else {
        compare(actual[key], expected[key]);
      }
    }
  } else if (typeof expected === 'number' && !Number.isInteger(expected)) {
    assert.ok(Math.abs(actual - expected) < 1e-10 && actual !== null);
  } else {
    assert.equal(actual, expected);
  }
}
for (const sample of fixtures) {
  test(sample.name, () => {
    const parsed = parseFasta(sample.text);
    compare(parsed, sample.parsed);
    compare(parsed.records.map(stats), sample.stats);
    compare(windows(parsed.records, 9, 5), sample.windows);
    parsed.records.forEach(({sequence}, index) => {
      const expected = sample.results[index];
      assert.equal(transcribe(sequence), expected.rna);
      for (const [table, frames] of Object.entries(expected.frames)) {
        compare(translate(sequence, Number(table)), frames);
      }
      for (const [options, hits] of expected.orfs) {
        compare(orfs(sequence, options), hits);
      }
      for (const [conditions, tm] of expected.tm) {
        compare({tm: melt(sequence, conditions)}, {tm});
      }
    });
  });
}
test('insulin default ORFs and invalid input', () => {
  const {sequence} = parseFasta(fixtures[0].text).records[0];
  const hits = orfs(sequence).map(hit => [hit.frame, hit.start, hit.end, hit.protein.length]);
  assert.ok(hits.some(hit => hit.join() === '+3,60,392,110'));
  assert.ok(hits.some(hit => hit.join() === '-3,206,355,49'));
  assert.throws(() => parseFasta('ACGT'));
  assert.throws(() => translate('AC-'));
  assert.throws(() => orfs('ATG', {table: 27}));
});
