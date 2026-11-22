const CACHE = 'dna-analyzer-v1';
const FILES = ['./', './index.html', './styles.css', './app.mjs', './core.mjs',
  './melting.mjs', './reference.mjs', './exports.mjs', './manifest.webmanifest',
  './icon.svg', './icon-192.png', './icon-512.png', './examples.fasta'];
const urls = new Set(FILES.map(file => new URL(file, self.registration.scope).href));

self.addEventListener('install', event => {
  event.waitUntil(caches.open(CACHE).then(cache => cache.addAll(
    FILES.map(file => new Request(file, {cache: 'reload'})),
  )));
});
self.addEventListener('activate', event => {
  event.waitUntil(caches.keys().then(keys => Promise.all(
    keys.filter(key => key.startsWith('dna-analyzer-') && key !== CACHE).map(key => caches.delete(key)),
  )).then(() => self.clients.claim()));
});
self.addEventListener('fetch', event => {
  if (event.request.method !== 'GET' || !urls.has(event.request.url)) return;
  event.respondWith(caches.open(CACHE).then(async cache =>
    (await cache.match(event.request)) ?? fetch(event.request),
  ));
});
