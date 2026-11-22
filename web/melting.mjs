import {complement, stacks} from './reference.mjs';

function saltCorrection(sequence, method, na, k, tris, magnesium, dntps) {
  if (!method) return 0;
  let monovalent = na + k + tris / 2;
  if (method !== 7 && dntps < magnesium) monovalent += 120 * Math.sqrt(magnesium - dntps);
  const mon = monovalent * 1e-3;
  if (method < 7 && !mon) throw new Error('La concentración total de sales debe ser positiva.');
  if (method === 1) return 16.6 * Math.log10(mon);
  if (method === 2) return 16.6 * Math.log10(mon / (1 + 0.7 * mon));
  if (method === 3) return 12.5 * Math.log10(mon);
  if (method === 4) return 11.7 * Math.log10(mon);
  if (method === 5) return 0.368 * (sequence.length - 1) * Math.log(mon);
  const gc = (sequence.match(/[GC]/g) ?? []).length / sequence.length;
  const logMon = Math.log(mon);
  let mg = magnesium * 1e-3;
  if (method === 7 && dntps > 0) {
    const association = 3e4;
    const term = association * (dntps * 1e-3 - mg) + 1;
    mg = (-term + Math.sqrt(term ** 2 + 4 * association * mg)) / (2 * association);
  }
  const ratio = Math.sqrt(mg) / mon;
  if (method === 6 || (mon > 0 && ratio < 0.22)) {
    return (4.29 * gc - 3.95) * 1e-5 * logMon + 9.40e-6 * logMon ** 2;
  }
  if (!mg) throw new Error('La concentración de magnesio libre debe ser positiva.');
  let a = 3.92;
  let d = 1.42;
  let g = 8.31;
  if (mon > 0 && ratio < 6) {
    a *= 0.843 - 0.352 * Math.sqrt(mon) * logMon;
    d *= 1.279 - 4.03e-3 * logMon - 8.03e-3 * logMon ** 2;
    g *= 0.486 - 0.258 * logMon + 5.25e-3 * logMon ** 3;
  }
  const logMg = Math.log(mg);
  return (a - 0.911 * logMg + gc * (6.26 + d * logMg)
    + (-48.2 + 52.5 * logMg + g * logMg ** 2) / (2 * (sequence.length - 1))) * 1e-5;
}

export function melt(sequence, {na_mM = 50, k_mM = 0, tris_mM = 0, mg_mM = 0, dntps_mM = 0,
  dnac1_nM = 25, dnac2_nM = 25, saltcorr = 5, selfcomp = false} = {}) {
  if (![na_mM, k_mM, tris_mM, mg_mM, dntps_mM, dnac1_nM, dnac2_nM].every(n => Number.isFinite(n) && n >= 0)) {
    throw new Error('Indica concentraciones finitas y no negativas.');
  }
  const strand = (dnac1_nM - (selfcomp ? 0 : dnac2_nM / 2)) * 1e-9;
  if (strand <= 0) throw new Error('La concentración efectiva de hebra debe ser positiva.');
  if (!Number.isInteger(saltcorr) || saltcorr < 0 || saltcorr > 7) throw new Error('Corrección salina entre 0 y 7.');
  const dna = sequence.replace(/\s/g, '').toUpperCase().replaceAll('U', 'T');
  if (!/^[ACGT]{2,}$/.test(dna)) return null;

  let [enthalpy, entropy] = stacks.init;
  for (const base of [dna[0], dna.at(-1)]) {
    const [h, s] = stacks['AT'.includes(base) ? 'init_A/T' : 'init_G/C'];
    enthalpy += h;
    entropy += s;
  }
  for (let i = 0; i + 1 < dna.length; i++) {
    const pair = dna.slice(i, i + 2) + '/' + complement[dna[i]] + complement[dna[i + 1]];
    const [h, s] = stacks[pair] ?? stacks[[...pair].reverse().join('')];
    enthalpy += h;
    entropy += s;
  }
  if (selfcomp) entropy += stacks.sym[1];
  const correction = saltCorrection(dna, saltcorr, na_mM, k_mM, tris_mM, mg_mM, dntps_mM);
  if (saltcorr === 5) entropy += correction;
  const kelvin = 1000 * enthalpy / (entropy + 1.987 * Math.log(strand));
  if (saltcorr >= 6) return 1 / (1 / kelvin + correction) - 273.15;
  return kelvin - 273.15 + (saltcorr < 5 ? correction : 0);
}
