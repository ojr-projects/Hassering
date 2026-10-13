# %% build_library.py: known bacterial motif library for the tendency-map viewer ("next nearest" matches)
"""
Sources (frequency matrices, rows A C G T per position):
  PRODORIC   every matrix of the PRODORIC API (https://www.prodoric.de/api/matrix); the matrix is rebuilt from its
             curated binding-site sequences (modal length, +0.1 pseudocount) when they exist, else from the published
             weights normalised per column
  RegulonDB  E. coli K-12 TF sites (RISet.gff3, riType tf-*), modal length, >= 5 sites
  classic    hand-entered bacterial elements from the literature (IUPAC; sources in CLASSIC below)
Output: tmap/motif_library.json  [{name, source, organism, consensus, n_sites, matrix}]
Usage: python tmap/build_library.py [regulondb_raw_dir]
RegulonDB sites are cut from the E. coli K-12 MG1655 genome (GCF_000005845.2); set ECOLI_GENOME to its FASTA.
"""
import collections
import json
import os
import sys
import time
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
API = 'https://www.prodoric.de/api'
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = str.maketrans('ACGT', 'TGCA')
GENOME = os.environ.get('ECOLI_GENOME', '/home/oliver/code/comsequent/data/genomes/Escherichia_coli_str_K-12_substr._MG1655/'
                        'ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna')
# classic elements: (name, IUPAC, organism/scope, source)
CLASSIC = [
    ('sigma70 -35 box', 'TTGACA', 'bacteria (sigma70)', 'Harley & Reynolds 1987'),
    ('sigma70 -10 box', 'TATAAT', 'bacteria (sigma70)', 'Harley & Reynolds 1987'),
    ('sigma70 -35..-10 (17 bp spacer)', 'TTGACANNNNNNNNNNNNNNNNNTATAAT', 'bacteria (sigma70)', 'Harley & Reynolds 1987'),
    ('extended -10 box', 'TGNTATAAT', 'bacteria (sigma70)', 'Mitchell et al. 2003'),
    ('sigma54 -24/-12', 'TGGCACGNNNNTTGCW', 'bacteria (sigma54)', 'Barrios et al. 1999'),
    ('sigma32 promoter -35/-10', 'TCTCNCCCTTGAANNNNNNNNNNNNNCCCCATNTA', 'E. coli (sigma32)', 'Nonaka et al. 2006'),
    ('Bacteroides promoter -33/-7', 'TTTGNNNNNNNNNNNNNNNNNNNNTANNTTTG', 'Bacteroides', 'Bayley et al. 2000; Mastropaolo 2009'),
    ('Shine-Dalgarno (RBS)', 'AGGAGG', 'bacteria', 'Shine & Dalgarno 1974'),
    ('Chi site (E. coli)', 'GCTGGTGG', 'Escherichia', 'Smith et al. 1981'),
    ('Chi site (Bacillus subtilis)', 'AGCGG', 'Bacillus', 'Chedin et al. 1998'),
    ('Chi site (Lactococcus lactis)', 'GCGCGTG', 'Lactococcus', 'Biswas et al. 1995'),
    ('Chi site (Staphylococcus aureus)', 'GAAGCGG', 'Staphylococcus', 'Halpern et al. 2007'),
    ('Chi site (Haemophilus influenzae)', 'GNTGGTGG', 'Haemophilus', 'Sourice et al. 1998'),
    ('DnaA box', 'TTWTNCACA', 'bacteria', 'Fuller et al. 1984; Schaper & Messer 1995'),
    ('HIP1 octamer', 'GCGATCGC', 'cyanobacteria', 'Robinson et al. 1995'),
    ('DNA uptake sequence (Neisseria DUS)', 'GCCGTCTGAA', 'Neisseria', 'Goodman & Scocca 1988'),
    ('DNA uptake sequence (Haemophilus USS)', 'AAGTGCGGT', 'Pasteurellaceae', 'Smith et al. 1995'),
    ('KOPS (FtsK orienting)', 'GGGNAGGG', 'E. coli', 'Bigot et al. 2005'),
    ('matS (MatP site)', 'GTGACRNYGTCAC', 'E. coli', 'Mercier et al. 2008'),
    ('dif site (XerCD)', 'GGTGCGCATAATGTATATTATGTTAAAT', 'E. coli', 'Blakely et al. 1993'),
    ('CRP site', 'TGTGANNNNNNTCACA', 'E. coli', 'Berg & von Hippel 1988'),
    ('LexA SOS box (E. coli)', 'CTGTNNNNNNNNACAG', 'E. coli', 'Walker 1984'),
    ('LexA SOS box (Gram-positive)', 'GAACNNNNGTTC', 'Bacillus / Firmicutes', 'Winterling et al. 1998'),
    ('Fur box', 'GATAATGATAATCATTATC', 'E. coli', 'de Lorenzo et al. 1987'),
    ('IHF site', 'WATCAANNNNTTR', 'E. coli', 'Goodrich et al. 1990'),
    ('Fis site', 'GNNNYRNNWNNNYRNNNC', 'E. coli', 'Shao et al. 2008'),
    ('ArcA site', 'GTTAATTAA', 'E. coli', 'Liu & De Wulf 2004'),
    ('FNR site', 'TTGATNNNNATCAA', 'E. coli', 'Spiro & Guest 1990'),
    ('CodY box', 'AATTTTCWGAAAATT', 'Firmicutes', 'den Hengst et al. 2005'),
    ('CcpA cre site', 'WTGNAANCGNWNNCW', 'Firmicutes', 'Weickert & Chambliss 1990'),
    ('rho-independent terminator U-tract', 'TTTTTTT', 'bacteria', 'd\'Aubenton Carafa et al. 1990'),
    ('REP element (E. coli PU)', 'GCCKGATGNCGRCGYNNNNNNRCGYCTTATCMGGCCTAC', 'E. coli', 'Gilson et al. 1984'),
    ('BIME box (E. coli)', 'GCCGGATG', 'E. coli', 'Bachellier et al. 1997'),
    ('T-box (tRNA antiterminator core)', 'AGGGUGGNACCGCG'.replace('U', 'T'), 'Firmicutes', 'Grundy & Henkin 1993'),
    ('PhoB pho box', 'CTGTCATANNNCTGTCAT', 'E. coli', 'Makino et al. 1988'),
    ('ArgR ARG box', 'TNTGAATNNAATTCANNA', 'E. coli', 'Maas 1994'),
    ('Lrp site', 'YAGHAWATTWTDCTR', 'E. coli', 'Cui et al. 1995'),
    ('TetR-family palindrome (generic)', 'TNNNNNNNNNNNNA', 'bacteria', 'placeholder generic dyad'),
]


def gff(path):
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(kv.split('=', 1) for kv in f[8].replace('; ', ';').split(';') if '=' in kv)
        yield f, a


def get(url, tries=4):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'tmap-library/1.0'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            time.sleep(2 * (t + 1))
    return None


def from_sites(seqs, pseudo=0.1):
    seqs = [s.upper() for s in seqs if s and set(s.upper()) <= set('ACGT')]
    if len(seqs) < 2:
        return None, 0
    L = collections.Counter(len(s) for s in seqs).most_common(1)[0][0]
    ss = [s for s in seqs if len(s) == L]
    C = np.array([[sum(s[i] == b for s in ss) for b in 'ACGT'] for i in range(L)], float) + pseudo
    return (C / C.sum(axis=1, keepdims=True)).round(4).tolist(), len(ss)


def from_iupac(s):
    M = []
    for ch in s.upper():
        allowed = IUPAC[ch]
        M.append([round(1 / len(allowed), 4) if b in allowed else 0.0 for b in 'ACGT'])
    return M


def consensus(M):
    out = ''
    for row in M:
        p = np.array(row)
        best = p.argmax()
        out += 'ACGT'[best] if p[best] >= 0.5 else 'acgt'[best] if p[best] >= 0.35 else 'n'
    return out


def main():
    lib = []
    lst = get(f'{API}/matrix') or []
    for k, m in enumerate(lst):
        if m.get('deprecated'):
            continue
        det = get(f'{API}/matrix/MX{m["mx"]:06d}') or {}
        M, n = from_sites([b.get('sequence', '') for b in det.get('bindingsites', [])])
        if M is None:
            W = np.array([m['pwm'][b] for b in 'ACGT'], float).T
            W = np.maximum(W, 0) + 0.05
            M, n = (W / W.sum(axis=1, keepdims=True)).round(4).tolist(), 0
        lib.append(dict(name=m['name'], source=f'PRODORIC MX{m["mx"]:06d}', organism=m['organism']['name'],
                        consensus=consensus(M), n_sites=n, matrix=M))
        if k % 50 == 0:
            print('prodoric', k, len(lst), flush=True)
    raw = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/regulondb/raw')
    if os.path.exists(os.path.join(raw, 'RISet.gff3')):
        g = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
        sites = collections.defaultdict(set)
        for f, a in gff(os.path.join(raw, 'RISet.gff3')):
            if a.get('riType', '').startswith('tf-'):
                sites[a['Name']].add((int(f[3]) - 1, int(f[4]), f[6]))
        for tf, ss in sorted(sites.items()):
            seqs = [g[a:b] if s == '+' else g[a:b].translate(COMP)[::-1] for a, b, s in ss]
            M, n = from_sites(seqs)
            if M is not None and n >= 5:
                lib.append(dict(name=tf, source='RegulonDB v12', organism='Escherichia coli K-12', consensus=consensus(M),
                                n_sites=n, matrix=M))
    for name, iu, org, src in CLASSIC:
        if 'placeholder' in src:
            continue
        M = from_iupac(iu)
        lib.append(dict(name=name, source=f'classic: {src}', organism=org, consensus=iu, n_sites=0, matrix=M))
    out = os.path.join(HERE, 'motif_library.json')
    json.dump(lib, open(out, 'w'), separators=(',', ':'))
    c = collections.Counter(e['source'].split(' ')[0].rstrip(':') for e in lib)
    print(f'{len(lib)} motifs -> {out} ({os.path.getsize(out) // 1024} KB): {dict(c)}')


if __name__ == '__main__':
    main()
