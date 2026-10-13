# %% bundles.py: cross-strain bundles of a family's occurrences (orthologous loci), pre-computed for the viewer
"""
Every broad-level occurrence of a family, in every strain, is assigned to a bundle in three priority stages:
  1. gene adjacency: the same flanking genes (Bakta gene name, else product; 'hypothetical protein' ignored) on the
     motif's upstream and downstream sides (oriented to the motif), at distances within DIST_TOL of the bundle's
  2. pattern context: the same surrounding structure, i.e. the set of (other family, offset/25 bp, relative strand)
     within +-WIN, Jaccard >= PAT_J with a bundle's representative
  3. sequence: MinHash (bottom NH of hashed 11-mers) of the oriented +-WIN sequence, Jaccard >= SEQ_J with a bundle's
     representative; otherwise the occurrence starts a new bundle
Bundles are ordered by size (strain count, then rows). Families with more than MAX_OCC occurrences are not bundled.
"""
import zlib

import numpy as np

WIN, DIST_TOL, PAT_J, SEQ_J, NH, K, MAX_OCC = 250, 80, 0.5, 0.2, 64, 11, 5000
RC = str.maketrans('ACGTN', 'TGCAN')
GENERIC = {'', 'hypothetical protein', 'putative protein', 'DUF domain-containing protein'}


def minhash(seq):
    s = seq.upper()
    if len(s) < K:
        return np.full(NH, np.iinfo(np.uint64).max, np.uint64)
    hs = np.unique(np.fromiter((zlib.crc32(s[i:i + K].encode()) for i in range(len(s) - K + 1)
                                if 'N' not in s[i:i + K]), np.uint64))          # stable across runs
    out = np.full(NH, np.iinfo(np.uint64).max, np.uint64)
    out[:min(NH, len(hs))] = hs[:NH]
    return out


def mh_sim(a, b):
    u = np.union1d(a, b)[:NH]
    return len(np.intersect1d(np.intersect1d(a, b), u)) / NH


def gene_label(f):
    """Bakta gene name when there is one (not the locus tag), else the product unless generic."""
    prod = (f[5] or '').strip()
    name = f[4] if f[4] and f[4] != f[6] else ''
    return name or (prod if prod not in GENERIC else '')


def flank(cds, lo, hi, strand):
    """(upstream_label, upstream_dist, downstream_label, downstream_dist), oriented to the motif."""
    import bisect
    starts = [c[0] for c in cds]
    i = bisect.bisect_left(starts, hi)
    right = cds[i] if i < len(cds) else None
    j = i - 1
    while j >= 0 and cds[j][1] > lo:                     # skip genes overlapping the motif
        j -= 1
    left = cds[j] if j >= 0 else None
    L = (gene_label(left), lo - left[1]) if left else ('', None)
    R = (gene_label(right), right[0] - hi) if right else ('', None)
    return (L[0], L[1], R[0], R[1]) if strand > 0 else (R[0], R[1], L[0], L[1])


def bundle_family(f, strains, genomes, nbr_index):
    """Returns (assignment {si: [bundle per occurrence]}, bundles [{rule, label, n, strains}])."""
    occ = [(si, k, o) for si, lst in f['occ'].items() for k, o in enumerate(lst) if o[2] >= f['broad_step']]
    if not occ or len(occ) > MAX_OCC:
        return None, None
    recs = []
    for si, k, (p, st, lv) in occ:
        si = int(si)
        lo = p if st > 0 else p - f['span'] + 1
        hi = lo + f['span']
        cds = strains[si]['cds_sorted']
        up, du, dn, dd = flank(cds, lo, hi, st)
        c = (lo + hi) // 2
        nb = set()
        for (a, b, s2, fid) in nbr_index[si](c - WIN, c + WIN):
            if fid == f['id'] and abs(a - lo) < 3:
                continue
            nb.add((fid, int(((a + b) // 2 - c) * st // 25), s2 * st))
        g = genomes[si]
        seq = g[max(0, c - WIN):c + WIN]
        if st < 0:
            seq = seq.translate(RC)[::-1]
        recs.append(dict(si=si, k=k, up=up, du=du, dn=dn, dd=dd, nb=nb, mh=minhash(seq)))
    bundles = []                                    # each: dict(rule, label, members, rep)
    # stage 1: gene adjacency
    by_key = {}
    for r in recs:
        if r['up'] or r['dn']:
            by_key.setdefault((r['up'], r['dn']), []).append(r)
    assigned = set()
    for key, rs in by_key.items():
        rs.sort(key=lambda r: (r['du'] or 0) + (r['dd'] or 0))
        groups = []
        for r in rs:
            for gr in groups:
                ref = gr[0]
                if (r['du'] is None or ref['du'] is None or abs(r['du'] - ref['du']) <= DIST_TOL) and \
                   (r['dd'] is None or ref['dd'] is None or abs(r['dd'] - ref['dd']) <= DIST_TOL):
                    gr.append(r)
                    break
            else:
                groups.append([r])
        for gr in groups:
            if len({r['si'] for r in gr}) >= 2:
                lab = f"{key[0] or '?'} ← motif → {key[1] or '?'}"
                bundles.append(dict(rule='gene adjacency', label=lab, members=gr, rep=gr[0]))
                assigned.update(id(r) for r in gr)
    # stage 2: pattern context; stage 3: sequence
    for r in recs:
        if id(r) in assigned:
            continue
        best, bs = None, 0.0
        if r['nb']:
            for b in bundles:
                nb2 = b['rep']['nb']
                if nb2:
                    j = len(r['nb'] & nb2) / len(r['nb'] | nb2)
                    if j >= PAT_J and j > bs:
                        best, bs = b, j
        if best is not None:
            best['members'].append(r)
            if best['rule'] == 'sequence':
                best['rule'] = 'pattern context'
            assigned.add(id(r))
            continue
        for b in bundles:
            s = mh_sim(r['mh'], b['rep']['mh'])
            if s >= SEQ_J and s > bs:
                best, bs = b, s
        if best is not None:
            best['members'].append(r)
        else:
            bundles.append(dict(rule='sequence', label='', members=[r], rep=r))
        assigned.add(id(r))
    for b in bundles:
        b['strains'] = len({r['si'] for r in b['members']})
        if len(b['members']) == 1:
            b['rule'] = 'singleton'
        if not b['label']:
            r = b['rep']
            b['label'] = f"{r['up'] or '?'} ← motif → {r['dn'] or '?'}"
    bundles.sort(key=lambda b: (-b['strains'], -len(b['members'])))
    assign = {str(si): [-1] * len(lst) for si, lst in f['occ'].items()}
    for bi, b in enumerate(bundles):
        for r in b['members']:
            assign[str(r['si'])][r['k']] = bi
    return assign, [dict(rule=b['rule'], label=b['label'], n=len(b['members']), strains=b['strains']) for b in bundles]
