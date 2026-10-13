# %% viewer_data.py: compact JSON for the circular tendency-map viewer (Revision 3)
"""
Reads <species>/spectra/<mode>.families.json (+ null_<mode>.zgap.npy / zballgap.npy: gap-stratified null for E-values) and each
strain's annotation (Bakta genome.gff3 when present, else Pyrodigal cds.bed). Writes <species>/spectra/<mode>.view.json:
  strains   [{acc, strain, length, records: [{name, offset, length}], feat: [[start, end, strand, type, label, product]]}]
  structs   [{id, families: [ids]}]
  families  [{id, struct, name, degenerate, dgroup, inverted, gap, span, z, e, z_ball, e_ball, n_seed, broad_step,
              chain: [[n, n_fixed]], cons_broad, cons_strict, lo_broad, lo_strict, ic_broad, mapback, bits (sum of
              information over the broad composite, bits), length, logo ([pA, pC, pG, pT] per column, composite +-4 bp,
              starting at offset logo_start), bg (domain base composition),
              occ: {strain_index: [[seed_start_fwd, strand, level], ...]}}]
Coordinates are forward-strand positions in the strain's concatenated record space (the circle). An occurrence's seed
start is the forward coordinate of the seed pair's first base on its own strand (strand -1: the pair reads leftward).
Usage: python tmap/viewer_data.py <species_dir> [mode=nc]
"""
import json
import os
import sys
import urllib.parse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats  # noqa: E402

LOGO_FLANK = 4
KEEP_TYPES = {'CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'regulatory_region', 'crispr-repeat', 'crispr-spacer',
              'repeat_region', 'oriC', 'oriT', 'ncRNA-region', 'sorf', 'region'}


def features(gdir, roff):
    gff = os.path.join(gdir, 'bakta', 'genome.gff3')
    out = []
    if os.path.exists(gff):
        for line in open(gff):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in KEEP_TYPES or f[0] not in roff:
                continue
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            lab = urllib.parse.unquote(a.get('gene') or a.get('Name') or a.get('locus_tag') or f[2])
            prod = urllib.parse.unquote(a.get('product', ''))
            out.append([roff[f[0]] + int(f[3]) - 1, roff[f[0]] + int(f[4]), 1 if f[6] == '+' else -1, f[2],
                        lab[:40], prod[:80], a.get('locus_tag', '')])
        return out, 'bakta'
    for line in open(os.path.join(gdir, 'cds.bed')):
        r, a, b, s = line.split()
        out.append([roff[r] + int(a), roff[r] + int(b), 1 if s == '+' else -1, 'CDS', '', '', ''])
    return out, 'pyrodigal'


def consensus(fm, step):
    span = fm['span']
    fx = {int(k): v for k, v in fm['chain'][step]['fixed'].items()}
    A, _, B = fm['name'].split('-')
    for i, ch in enumerate(A):
        fx[i] = ch
    for i, ch in enumerate(B):
        fx[span - len(B) + i] = ch
    lo, hi = min(fx), max(fx)
    return lo, ''.join(fx.get(o, '.') for o in range(lo, hi + 1))


def main():
    d = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'nc'
    fam = json.load(open(os.path.join(d, 'spectra', f'{mode}.families.json')))
    nz = os.path.join(d, 'spectra', f'null_{mode}.zgap.npy')           # gap-stratified null (stats.null_by_gap)
    nzb = os.path.join(d, 'spectra', f'null_{mode}.zballgap.npy')
    null = np.load(nz) if os.path.exists(nz) else None
    nullb = np.load(nzb) if os.path.exists(nzb) else None
    meta = {}
    pt = os.path.join(d, 'panel.tsv')
    if os.path.exists(pt):
        rows = open(pt).read().splitlines()
        hdr = rows[0].split('\t')
        for r in rows[1:]:
            v = dict(zip(hdr, r.split('\t')))
            meta[v['acc']] = v
    strains, src = [], set()
    for p in fam['strains']:
        gdir = os.path.dirname(p)
        acc = os.path.basename(gdir)
        recs, off, name, n = [], 0, None, 0
        for line in open(os.path.join(gdir, 'genome.fna')):
            if line.startswith('>'):
                if name is not None:
                    recs.append(dict(name=name, offset=off, length=n))
                    off += n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
        recs.append(dict(name=name, offset=off, length=n))
        roff = {r['name']: r['offset'] for r in recs}
        feat, s = features(gdir, roff)
        src.add(s)
        strains.append(dict(acc=acc, strain=meta.get(acc, {}).get('strain', ''), length=off + n, records=recs, feat=feat))
    rlen = {(i, r['name']): r['length'] for i, s in enumerate(strains) for r in s['records']}
    roff = {(i, r['name']): r['offset'] for i, s in enumerate(strains) for r in s['records']}
    gc = 0.0                                                     # background base composition of the active domain
    cnt = np.zeros(4)
    for p in fam['strains']:
        x = open(p).read()
        for k, ch in enumerate('ACGT'):
            cnt[k] += x.count(ch)
    bg = [round(float(v), 4) for v in cnt / cnt.sum()]
    Ppath = os.path.join(d, 'spectra', f'{mode}.pairs.u32')
    ZB = stats.ball_z(np.fromfile(Ppath, np.uint32).reshape(97, 256, 256))
    fams = []
    for fm in fam['families']:
        g0, a0, b0 = fm['ball_cells']
        fm['z_ball'] = round(float(ZB[g0, a0, b0]), 1)
        lob, cb = consensus(fm, fm['broad_step'])
        los, cs = consensus(fm, len(fm['chain']) - 1)
        L = fam['flank']
        prof = fm['chain'][fm['broad_step']]['profile']
        ic = [round(c[4], 2) for c in prof][max(0, L + lob):L + lob + len(cb)]
        # sequence logo: base frequencies of the aligned broad-level occurrences over the composite +- LOGO_FLANK
        l0, l1 = max(0, L + lob - LOGO_FLANK), min(len(prof), L + lob + len(cb) + LOGO_FLANK)
        logo = [c[:4] for c in prof[l0:l1]]
        bits = round(float(sum(ic)), 1)
        occ, seen = {}, {}
        for si, name, st, i, lv in sorted(fm['occ'], key=lambda o: o[2] != '+'):      # + strand first
            n = rlen[(si, name)]
            pos = i if st == '+' else n - 1 - i                  # forward coordinate of the seed's first base
            p0 = roff[(si, name)] + pos
            lo = p0 if st == '+' else p0 - fm['span'] + 1
            # an inverted (palindromic) pair matches the same physical site on both strands: keep one copy per footprint
            key = (si, lo)
            if key in seen:
                seen[key][2] = max(seen[key][2], lv)
                continue
            rec = [p0, 1 if st == '+' else -1, lv]
            seen[key] = rec
            occ.setdefault(si, []).append(rec)
        for lst in occ.values():
            lst.sort(key=lambda o: o[0])
        n_seed_sites = sum(len(v) for v in occ.values())
        fams.append(dict(
            id=fm['id'], struct=fm['structure'], name=fm['name'], degenerate=fm['degenerate'], dgroup=fm['degenerate_group'],
            inverted=fm['inverted'], gap=fm['gap'], span=fm['span'], z=fm['z'], z_ball=fm['z_ball'],
            count=fm['count'], baseline=round(fm['baseline'], 1), count_ball=fm['count_ball'],
            e=None if null is None else stats.e_value_gap(fm['z'], fm['gap'], null),
            e_ball=None if nullb is None else stats.e_value_gap(fm['z_ball'], fm['gap'], nullb),
            n_seed=n_seed_sites, n_seed_strand_hits=fm['n_seed'], broad_step=fm['broad_step'], chain=[[c['n'], len(c['fixed'])] for c in fm['chain']],
            cons_broad=cb, lo_broad=lob, cons_strict=cs, lo_strict=los, ic_broad=ic, mapback=fm['mapback'], occ=occ,
            bits=bits, length=len(cb), logo=logo, logo_start=l0 - L, bg=bg))
    # cross-strain bundles (tmap/bundles.py) for every family a default viewer could show (E <= 1e-2)
    import bisect
    import bundles as bmod
    genomes = []
    for st in strains:
        p = os.path.join(os.path.dirname(fam['strains'][strains.index(st)]), 'genome.fna')
        genomes.append(''.join(l.strip().upper() for l in open(p) if not l.startswith('>')))
        st['cds_sorted'] = sorted((f for f in st['feat'] if f[3] == 'CDS'), key=lambda f: f[0])
    allocc = [[] for _ in strains]
    for f in fams:
        for si, lst in f['occ'].items():
            for p0, s0, lv in lst:
                if lv >= f['broad_step']:
                    lo = p0 if s0 > 0 else p0 - f['span'] + 1
                    allocc[int(si)].append((lo, lo + f['span'], s0, f['id']))
    for a_ in allocc:
        a_.sort()
    starts = [[x[0] for x in a_] for a_ in allocc]

    def make_query(si):
        def q(lo, hi):
            i = bisect.bisect_left(starts[si], lo - 200)
            out = []
            while i < len(allocc[si]) and allocc[si][i][0] <= hi:
                if allocc[si][i][1] >= lo:
                    out.append(allocc[si][i])
                i += 1
            return out
        return q
    nbr = [make_query(si) for si in range(len(strains))]
    nb_done = 0
    for f in fams:
        if f['e'] is not None and f['e'] <= 1e-2:
            f['bundle'], f['bundles'] = bmod.bundle_family(f, strains, genomes, nbr)
            nb_done += f['bundle'] is not None
    for st in strains:
        st.pop('cds_sorted')
    print(f'bundled {nb_done} families', flush=True)
    structs = [dict(id=k, families=g) for k, g in enumerate(fam['structures'])]
    out = dict(species=os.path.basename(d).replace('_', ' '), mode=mode, annotation=sorted(src), strains=strains,
               structs=structs, families=fams, null_cells=None if null is None else int(null.size))
    p = os.path.join(d, 'spectra', f'{mode}.view.json')
    json.dump(out, open(p, 'w'), separators=(',', ':'))
    print(p, os.path.getsize(p) // 1024, 'KB')


if __name__ == '__main__':
    main()
