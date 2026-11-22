import {dna} from './core.mjs';

const MAX_URL_LENGTH = 2000;

export function blastUrl(sequence, program = 'blastn') {
  const query = program === 'blastn' ? dna(sequence) : sequence;
  if (!query) throw new Error('Selecciona un registro con secuencia.');
  const url = new URL('https://blast.ncbi.nlm.nih.gov/Blast.cgi');
  url.searchParams.set('PAGE_TYPE', 'BlastSearch');
  url.searchParams.set('PROGRAM', program);
  if (program === 'blastn') url.searchParams.set('MEGABLAST', 'on');
  url.searchParams.set('QUERY', query);
  if (url.href.length > MAX_URL_LENGTH) url.searchParams.delete('QUERY');
  return url.href;
}
