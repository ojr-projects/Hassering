# %% serve.py: local server for the tendency-map viewer (no artifact size limit; sequence on demand)
"""
GET  /                          the viewer (tmap/viewer_template.html, data fetched from the API)
GET  /api/species               [{name, modes}] for every panel with built view files
GET  /api/library               the known-motif library (tmap/motif_library.json)
GET  /api/view?species=NAME     {species, strains, views: {mode: view}}: the full data, no sampling cap
POST /api/seqs                  body [[species, strain_index, start, end], ...] -> ["ACGT...", ...] forward-strand slices
                                of the strain's concatenated records (viewer coordinates)
Usage: python tmap/serve.py [port=8765] [panels_root=~/data/tmap/panels]
Open http://localhost:<port>/ (WSL2 forwards localhost to Windows browsers).
"""
import json
import os
import sys
import urllib.parse
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
ROOT = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else '~/data/tmap/panels')
MODES = ('nc', 'all')


def species_list():
    out = []
    for sp in sorted(os.listdir(ROOT)):
        d = os.path.join(ROOT, sp, 'spectra')
        m = [x for x in MODES if os.path.exists(os.path.join(d, f'{x}.view.json'))]
        if m:
            out.append(dict(name=sp, label=sp.replace('_', ' '), modes=m))
    return out


@lru_cache(maxsize=4)
def view(sp, mtime_key):
    views = {}
    for m in MODES:
        p = os.path.join(ROOT, sp, 'spectra', f'{m}.view.json')
        if os.path.exists(p):
            views[m] = json.load(open(p))
    strains = next(iter(views.values()))['strains']
    for v in views.values():
        v.pop('strains')
        for f in v['families']:
            lv = {'broad': f['broad_step'], 'strict': len(f['chain']) - 1}
            f['true_counts'] = {k: [sum(1 for o in f['occ'].get(str(i), []) if o[2] >= m) for i in range(len(strains))]
                                for k, m in lv.items()}
            f['sampled'] = 1
    return json.dumps(dict(species=sp.replace('_', ' '), strains=strains, views=views), separators=(',', ':')).encode()


@lru_cache(maxsize=64)
def genome(sp, si):
    fam = json.load(open(os.path.join(ROOT, sp, 'spectra', 'nc.families.json')))
    p = os.path.join(os.path.dirname(fam['strains'][si]), 'genome.fna')
    return ''.join(l.strip().upper() for l in open(p) if not l.startswith('>'))


class H(BaseHTTPRequestHandler):
    def _send(self, body, ctype='application/json', code=200):
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *a):
        pass

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if u.path in ('/', '/index.html'):
            self._send(open(os.path.join(HERE, 'viewer_template.html'), 'rb').read(), 'text/html; charset=utf-8')
        elif u.path == '/api/library':
            lib = os.path.join(HERE, 'motif_library.json')
            self._send(open(lib, 'rb').read() if os.path.exists(lib) else b'[]')
        elif u.path == '/api/species':
            self._send(json.dumps(species_list()).encode())
        elif u.path == '/api/view':
            sp = q.get('species', [''])[0]
            d = os.path.join(ROOT, sp, 'spectra')
            if not os.path.isdir(d):
                return self._send(b'{"error":"unknown species"}', code=404)
            key = tuple(os.path.getmtime(os.path.join(d, f'{m}.view.json')) for m in MODES
                        if os.path.exists(os.path.join(d, f'{m}.view.json')))
            self._send(view(sp, key))
        else:
            self._send(b'{"error":"not found"}', code=404)

    def do_POST(self):
        if self.path != '/api/seqs':
            return self._send(b'{"error":"not found"}', code=404)
        req = json.loads(self.rfile.read(int(self.headers.get('Content-Length', 0))) or b'[]')
        out = []
        for sp, si, a, b in req[:2000]:
            g = genome(sp, int(si))
            out.append(g[max(0, int(a)):max(0, int(b))])
        self._send(json.dumps(out).encode())


if __name__ == '__main__':
    print(f'tendency map viewer on http://localhost:{PORT}/  (panels: {ROOT})', flush=True)
    ThreadingHTTPServer(('127.0.0.1', PORT), H).serve_forever()
