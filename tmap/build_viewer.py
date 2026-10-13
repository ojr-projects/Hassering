# %% build_viewer.py: inject a species' view JSON (every available mode) into the viewer template -> standalone HTML
"""Usage: python tmap/build_viewer.py <species_dir> <out.html> [modes=nc,all]"""
import json
import os
import sys

d, out = sys.argv[1], sys.argv[2]
modes = (sys.argv[3] if len(sys.argv) > 3 else 'nc,all').split(',')
views = {}
for m in modes:
    p = os.path.join(d, 'spectra', f'{m}.view.json')
    if os.path.exists(p):
        views[m] = json.load(open(p))
CAP = 4000          # occurrences drawn per family; exact per-strain counts are kept for every family
for v in views.values():
    for f in v['families']:
        lv = {'broad': f['broad_step'], 'strict': len(f['chain']) - 1}
        f['true_counts'] = {k: [sum(1 for o in f['occ'].get(str(i), []) if o[2] >= m) for i in range(len(v['strains']))]
                            for k, m in lv.items()}
        tot = sum(len(x) for x in f['occ'].values())
        k = -(-tot // CAP)
        f['sampled'] = k
        if k > 1:
            f['occ'] = {si: o[::k] for si, o in f['occ'].items()}
species = next(iter(views.values()))['species']
strains = next(iter(views.values()))['strains']           # identical in every mode: stored once
for v in views.values():
    v.pop('strains')
lib = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'motif_library.json')))
data = json.dumps(dict(species=species, strains=strains, views=views, library=lib), separators=(',', ':')).replace('</', '<\\/')
tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viewer_template.html')).read()
open(out, 'w').write(tpl.replace('__DATA__', data))
print(out, os.path.getsize(out) // 1024, 'KB', 'modes', list(views))
