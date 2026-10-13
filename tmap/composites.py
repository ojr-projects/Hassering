# %% composites.py: exact composite families from the offset-pair spectrum (Revisions 2A, 3)
"""
1. Seed pairs: the TOP canonical (a, g, b) cells of the spacing spectrum by spacing z (stats.spacing_z). The pair and its
   reverse-complement view are one family; the compact name is a-g-b (e.g. TTAG-96-CTAA).
2. Exact occurrences of the seed: core a at i and core b at i + c + g, footprint N-free in the active domain, both
   strands (oriented to the pair).
3. Growth (exact Apriori chain): align occurrences with L flanking bases; repeatedly add the position-base shared by
   the most current occurrences (>= KEEP of them). Each step is an exact composite whose occurrence set is a subset of
   its parent's. Every occurrence records the number of steps it survives (lv), so any level of the chain can be shown:
     broad level  = the deepest step keeping >= BROAD of the seed occurrences (default display)
     strict level = the end of the chain
4. Degenerate form: the Hamming-1 ball of both cores at the same gap (stats.ball_z, stats.iupac_ball): ball count,
   ball z and the IUPAC seed form (e.g. TTRG-96-CYAA). Families whose balls share cells form a degenerate group.
5. Structures: families whose broad-level occurrence loci overlap (Jaccard >= JAC) are one structure.
6. Map-back: the seed pair and the broad composite are counted exactly in every strain's whole genome, split by
   whether the footprint lies in CDS (cds.bed) or not.
Output: <species>/spectra/<mode>.families.json
Usage: python tmap/composites.py <species_dir> [mode=nc] [TOP=300]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats  # noqa: E402

C, GMAX, L, KEEP, BROAD, JOIN, JAC = 4, 96, 40, 0.5, 0.9, 10, 0.3
ENC = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'):
    ENC[ch] = i
    ENC[ch + 32] = i
RC = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(RC)[::-1]


def read_fasta(p):
    recs, name, buf = [], None, []
    for line in open(p):
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def core_codes(x):
    n = len(x) - C + 1
    if n <= 0:
        return np.zeros(0, np.int64), np.zeros(0, bool)
    code = np.zeros(n, np.int64)
    bad = np.zeros(n, bool)
    for j in range(C):
        col = x[j:j + n]
        bad |= col > 3
        code = code * 4 + np.where(col > 3, 0, col)
    return code, ~bad


def next_n(x):
    n = len(x)
    isn = np.flatnonzero(x > 3)
    if len(isn) == 0:
        return np.full(n, n)
    idx = np.searchsorted(isn, np.arange(n))
    return np.where(idx < len(isn), isn[np.minimum(idx, len(isn) - 1)], n)


def load_strands(paths):
    """[(strain, record, strand, seq_str, codes, core, okcore, nextN)] for both orientations."""
    out = []
    for si, p in enumerate(paths):
        for name, s in read_fasta(p):
            for st, x in (('+', s), ('-', rc(s))):
                e = ENC[np.frombuffer(x.encode(), np.uint8)]
                code, ok = core_codes(e)
                out.append((si, name, st, x, e, code, ok, next_n(e)))
    return out


def seed_hits(S, a, g, b):
    occ = []
    for si, name, st, x, e, code, ok, nn in S:
        n = len(code)
        if n == 0:
            continue
        j = np.arange(n) + C + g
        m = (j < n) & ok
        jj = np.where(m, j, 0)
        m &= ok[jj] & (code == a) & (code[jj] == b) & (j + C <= nn[:n])
        occ += [(si, name, st, int(i)) for i in np.flatnonzero(m)]
    return occ


def grow(index, occ, span, bg):
    """Exact Apriori chain. Returns chain steps and, per seed occurrence, the number of steps it survives."""
    rows = []
    for si, name, st, i in occ:
        x = index[(si, name, st)]
        w = x[max(0, i - L):i + span + L]
        pad = L - (i - max(0, i - L))
        rows.append('N' * pad + w + 'N' * (span + 2 * L - pad - len(w)))
    M = np.array([np.frombuffer(r.encode(), np.uint8) for r in rows]) if rows else np.zeros((0, span + 2 * L), np.uint8)
    alive = np.ones(len(rows), bool)
    lv = np.zeros(len(rows), np.int32)
    fixed, chain = {}, []
    while True:
        cur = M[alive]
        prof = []
        for col in range(M.shape[1]):
            f = np.array([(cur[:, col] == ord(ch)).sum() for ch in 'ACGT'], float)
            n = f.sum()
            p = f / n if n else np.full(4, .25)
            ic = float(np.sum(p[p > 0] * np.log2(p[p > 0] / bg[p > 0]))) if n else 0.0
            prof.append([round(float(v), 3) for v in p] + [round(ic, 3)])
        chain.append(dict(fixed={str(k): v for k, v in sorted(fixed.items())}, n=int(alive.sum()), profile=prof))
        best = None
        for col in range(M.shape[1]):
            o = col - L
            if o in fixed or 0 <= o < C or span - C <= o < span:
                continue
            vals, cnt = np.unique(cur[:, col], return_counts=True)
            for v, k in zip(vals, cnt):
                if v != ord('N') and k >= KEEP * len(cur) and (best is None or k > best[2]):
                    best = (o, chr(v), k)
        if best is None or len(fixed) >= 30:
            return chain, lv
        o, v, _ = best
        fixed[o] = v
        keep = alive & (M[:, o + L] == ord(v))
        lv[keep] += 1
        alive = keep


def composite_hits(Sall, fixed, span):
    """Exact occurrences of a composite {offset: base} (offsets relative to the seed start) in whole-genome strands."""
    items = sorted(fixed.items())
    lo = min(o for o, _ in items)
    hi = max(o for o, _ in items)
    out = []
    for si, name, st, x, e, code, ok, nn in Sall:
        n = len(e)
        idx = np.arange(max(0, -lo), n - hi)
        if len(idx) == 0:
            continue
        m = np.ones(len(idx), bool)
        for o, ch in items:
            m &= e[idx + o] == 'ACGT'.index(ch)
            if not m.any():
                break
        out += [(si, name, st, int(i)) for i in idx[m]]
    return out


def main():
    d = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'nc'
    top = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    paths = [l.strip() for l in open(os.path.join(d, f'genomes_{mode}.txt')) if l.strip()]
    P = np.fromfile(os.path.join(d, 'spectra', f'{mode}.pairs.u32'), np.uint32).reshape(GMAX + 1, 4 ** C, 4 ** C)
    Pf = P.astype(np.float64)
    base = stats.spacing_base(Pf)
    Z = (Pf - base) / np.sqrt(base + 1.0)
    ZB = stats.ball_z(Pf)
    B = stats.ball_counts(Pf)
    S = load_strands(paths)
    index = {(si, name, st): x for si, name, st, x, *_ in S}
    allx = ''.join(x for _, _, st, x, *_ in S if st == '+')
    bg = np.array([allx.count(ch) for ch in 'ACGT'], float)
    bg /= bg.sum()
    # whole-genome strands and CDS masks for map-back
    gpaths = [os.path.join(os.path.dirname(p), 'genome.fna') for p in paths]
    Sall = load_strands(gpaths)
    cdsmask = {}
    for si, p in enumerate(gpaths):
        lens = {n: len(s) for n, s in read_fasta(p)}
        for n, ln in lens.items():
            cdsmask[(si, n)] = np.zeros(ln, bool)
        for line in open(os.path.join(os.path.dirname(p), 'cds.bed')):
            r, a0, b0, _ = line.split()
            cdsmask[(si, r)][int(a0):int(b0)] = True
    rlen = {(si, n): len(m) for (si, n), m in cdsmask.items()}

    def in_cds(si, name, st, i, width):
        n = rlen[(si, name)]
        a0 = i if st == '+' else n - i - width
        return bool(cdsmask[(si, name)][max(0, a0):a0 + width].any())

    seen, fams = set(), []
    for f in np.argsort(Z, axis=None)[::-1]:
        g, a, b = (int(v) for v in np.unravel_index(f, Z.shape))
        A, Bc = stats.dec(a), stats.dec(b)
        key = min((A, Bc, g), (rc(Bc), rc(A), g))
        if key in seen:
            continue
        seen.add(key)
        span = 2 * C + g
        occ = seed_hits(S, a, g, b)
        chain, lv = grow(index, occ, span, bg)
        n0 = chain[0]['n']
        broad = max(s for s, c in enumerate(chain) if c['n'] >= BROAD * n0)
        ia, ib = stats.iupac_ball(Pf, base[g], g, a, b)
        # map-back of the seed pair and the broad composite into whole genomes (coding vs non-coding footprint)
        fx_b = {int(k): v for k, v in chain[broad]['fixed'].items()}
        for i, ch in enumerate(A):
            fx_b[i] = ch
        for i, ch in enumerate(Bc):
            fx_b[span - C + i] = ch
        lo, hi = min(fx_b), max(fx_b)
        mb = {}
        seed = seed_hits(Sall, a, g, b)            # every broad occurrence contains the seed pair
        extra = {o: ch for o, ch in fx_b.items() if not (0 <= o < C or span - C <= o < span)}
        idx_all = {(si, name, st): e for si, name, st, x, e, *_ in Sall}
        broad_hits = []
        for si, name, st, i in seed:
            e = idx_all[(si, name, st)]
            if all(0 <= i + o < len(e) and e[i + o] == 'ACGT'.index(ch) for o, ch in extra.items()):
                broad_hits.append((si, name, st, i))
        for lab, hits, off, w in (('seed', seed, 0, span), ('broad', broad_hits, lo, hi - lo + 1)):
            cod = sum(in_cds(si, name, st, i + off, w) for si, name, st, i in hits)
            mb[lab] = dict(total=len(hits), coding=cod, noncoding=len(hits) - cod)
        fams.append(dict(id=len(fams), name=f'{A}-{g}-{Bc}', degenerate=f'{ia}-{g}-{ib}', inverted=Bc == rc(A), gap=g,
                         z=round(float(Z[g, a, b]), 1), count=int(P[g, a, b]), baseline=float(base[g, a, b]),
                         z_ball=round(float(ZB[g, a, b]), 1), count_ball=int(B[g, a, b]), ball_cells=[g, a, b],
                         span=span, n_seed=n0, broad_step=broad, chain=chain,
                         occ=[[o[0], o[1], o[2], o[3], int(l)] for o, l in zip(occ, lv)], mapback=mb))
        if len(fams) >= top:
            break
        if len(fams) % 50 == 0:
            print(len(fams), 'families', flush=True)
    # degenerate groups: families whose Hamming balls share a cell (same gap, cores within 2 mismatches pairwise)
    def ham(x, y):
        return sum(((x >> (2 * j)) & 3) != ((y >> (2 * j)) & 3) for j in range(C))
    par = list(range(len(fams)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for i in range(len(fams)):
        gi, ai, bi = fams[i]['ball_cells']
        for j in range(i + 1, len(fams)):
            gj, aj, bj = fams[j]['ball_cells']
            if gi == gj and ham(ai, aj) <= 2 and ham(bi, bj) <= 2 and (ham(ai, aj) + ham(bi, bj)) <= 2:
                par[find(i)] = find(j)
    for i, fm in enumerate(fams):
        fm['degenerate_group'] = find(i)
    # structures from broad-level loci
    def loci(fm):
        out = set()
        for si, name, st, i, l in fm['occ']:
            if l < fm['broad_step']:
                continue
            n = rlen[(si, name)]
            start = i if st == '+' else n - i - fm['span']
            out.add((si, name, start // JOIN))
        return out
    LS = [loci(fm) for fm in fams]
    par2 = list(range(len(fams)))

    def find2(x):
        while par2[x] != x:
            par2[x] = par2[par2[x]]
            x = par2[x]
        return x
    for i in range(len(fams)):
        for j in range(i + 1, len(fams)):
            inter = len(LS[i] & LS[j])
            if inter and inter / len(LS[i] | LS[j]) >= JAC:
                par2[find2(i)] = find2(j)
    groups = {}
    for i in range(len(fams)):
        groups.setdefault(find2(i), []).append(i)
    structs = sorted(groups.values(), key=lambda g: -max(len(LS[i]) for i in g))
    for k, g in enumerate(structs):
        for i in g:
            fams[i]['structure'] = k
    out = dict(mode=mode, n_strains=len(paths), strains=paths, core=C, gmax=GMAX, flank=L, keep=KEEP, broad=BROAD,
               families=fams, structures=[sorted(g, key=lambda i: -len(LS[i])) for g in structs])
    json.dump(out, open(os.path.join(d, 'spectra', f'{mode}.families.json'), 'w'))
    print(f'{len(fams)} families in {len(structs)} structures -> {mode}.families.json')


if __name__ == '__main__':
    main()
