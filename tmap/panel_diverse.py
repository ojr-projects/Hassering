# %% panel_diverse.py: diversity-selected strain panels (MinHash/Mash distance, farthest-point selection)
"""
For a species with more complete non-MAG genomes than the panel size: download every complete genome's FASTA into
<root>/<species>/_pool/, sketch each one (bottom-S MinHash of canonical 21-mers), and choose N by farthest-point
selection. The selection starts from the medoid and skips genomes closer than MIN_D to an already-chosen one
(clones). The panel directory then holds exactly the chosen genomes; earlier picks not chosen move to _pool.
panel.tsv is rewritten with each genome's distance to its nearest other chosen genome.
Usage: python tmap/panel_diverse.py <root> <species> [N=15]
"""
import json
import os
import shutil
import sys
import urllib.parse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_fetch import API, get  # noqa: E402

K, S, MIN_D = 21, 2000, 0.0005
ENC = np.full(256, 255, np.uint8)
for i, c in enumerate(b'ACGT'):
    ENC[c] = i


def read_fna(p):
    seqs, cur = [], []
    for line in open(p, 'rb'):
        if line.startswith(b'>'):
            if cur:
                seqs.append(b''.join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append(b''.join(cur))
    return seqs


def sketch(p):
    hs = []
    mask = (1 << (2 * K)) - 1
    for s in read_fna(p):
        x = ENC[np.frombuffer(s, np.uint8)]
        n = len(x) - K + 1
        if n <= 0:
            continue
        bad = np.convolve((x > 3).astype(np.int32), np.ones(K, np.int32), 'valid') > 0
        xx = np.where(x > 3, 0, x).astype(np.uint64)
        f = np.zeros(n, np.uint64)
        r = np.zeros(n, np.uint64)
        for j in range(K):
            f = (f << np.uint64(2)) | xx[j:j + n]
            r = r | ((np.uint64(3) - xx[j:j + n]) << np.uint64(2 * j))
        c = np.minimum(f & np.uint64(mask), r & np.uint64(mask))[~bad]
        h = c * np.uint64(0x9E3779B97F4A7C15)
        h ^= h >> np.uint64(29)
        hs.append(np.unique(h))
    h = np.unique(np.concatenate(hs))
    return h[:S] if len(h) > S else h


def mash(a, b):
    u = np.union1d(a, b)[:S]
    j = len(np.intersect1d(np.intersect1d(a, b), u)) / len(u)
    return 1.0 if j == 0 else -np.log(2 * j / (1 + j)) / K


def main():
    root, sp = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    d = os.path.join(root, sp.replace(' ', '_'))
    pool = os.path.join(d, '_pool')
    os.makedirs(pool, exist_ok=True)
    q = urllib.parse.quote(sp)
    reps, token = [], None
    while True:
        r = get(f'{API}/genome/taxon/{q}/dataset_report?filters.is_metagenome_derived=METAGENOME_DERIVED_EXCLUDE'
                f'&filters.exclude_atypical=true&filters.assembly_level=complete_genome&page_size=1000'
                + (f'&page_token={token}' if token else ''))
        reps += r.get('reports', [])
        token = r.get('next_page_token')
        if not token:
            break
    accs, seen = [], set()
    for g in reps:
        key = g['accession'].split('_')[1].split('.')[0]
        if key not in seen:
            seen.add(key)
            accs.append((g['accession'], g.get('organism', {}).get('infraspecific_names', {}).get('strain', '')))
    import io
    import zipfile
    paths = {}
    for a, _ in accs:
        fp = os.path.join(d, a, 'genome.fna') if os.path.exists(os.path.join(d, a, 'genome.fna')) else os.path.join(pool, a + '.fna')
        if not os.path.exists(fp):
            z = get(f'{API}/genome/accession/{a}/download?include_annotation_type=GENOME_FASTA', binary=True)
            with zipfile.ZipFile(io.BytesIO(z)) as zf:
                m = [x for x in zf.namelist() if x.endswith('.fna')][0]
                open(fp, 'wb').write(zf.read(m))
        paths[a] = fp
    names = [a for a, _ in accs]
    sk = {a: sketch(paths[a]) for a in names}
    m = len(names)
    D = np.zeros((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            D[i, j] = D[j, i] = mash(sk[names[i]], sk[names[j]])
    chosen = [int(np.argmin(D.sum(axis=1)))]
    while len(chosen) < min(n, m):
        dmin = D[:, chosen].min(axis=1)
        dmin[chosen] = -1
        k = int(np.argmax(dmin))
        if dmin[k] < MIN_D:
            break
        chosen.append(k)
    sel = {names[i] for i in chosen}
    for a in os.listdir(d):                       # move unchosen earlier picks out of the panel
        if a.startswith('GC') and a not in sel and os.path.isdir(os.path.join(d, a)):
            shutil.move(os.path.join(d, a), os.path.join(pool, a))
    strain = dict(accs)
    with open(os.path.join(d, 'panel.tsv'), 'w') as fh:
        fh.write('acc\tstrain\tlevel\tnearest_chosen_mash\n')
        for i in chosen:
            a = names[i]
            gd = os.path.join(d, a)
            os.makedirs(gd, exist_ok=True)
            if not os.path.exists(os.path.join(gd, 'genome.fna')):
                src = paths[a] if os.path.exists(paths[a]) else os.path.join(pool, a, 'genome.fna')
                shutil.copy(src, os.path.join(gd, 'genome.fna'))
            near = min(D[i, j] for j in chosen if j != i) if len(chosen) > 1 else 0
            fh.write(f'{a}\t{strain[a]}\tComplete Genome\t{near:.4f}\n')
    json.dump(dict(species=sp, n_complete=m, chosen=[names[i] for i in chosen],
                   mash_range=[float(D[np.ix_(chosen, chosen)][np.triu_indices(len(chosen), 1)].min()),
                               float(D[np.ix_(chosen, chosen)].max())],
                   species_max=float(D.max())), open(os.path.join(d, 'diversity.json'), 'w'), indent=1)
    sub = D[np.ix_(chosen, chosen)][np.triu_indices(len(chosen), 1)]
    print(f'{sp}: {m} complete genomes; chosen {len(chosen)}; pairwise Mash among chosen {sub.min():.4f}-{sub.max():.4f} '
          f'(species max {D.max():.4f})', flush=True)


if __name__ == '__main__':
    main()
