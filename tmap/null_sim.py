# %% null_sim.py: Markov-simulated null panel for tendency-map E-values
"""
Fits an order-K Markov chain (K = 5) to the active sequence domain of a species panel (both strands, pooled), then
writes a synthetic twin of every strain: each run of non-N sequence is regenerated from the chain at the same length,
and N positions (masked coding sequence, record ends) stay where they are. The synthetic panel has the panel's
composition and short-range context but no longer-range structure: no repeats, no fixed-spacing pairs.
Then runs `tmap spectrum` on it and saves the null spacing-peak statistics (pairs z, and ball z) for E-values:
  <species>/spectra/null_<mode>.zsorted.npy        sorted z of every (g, a, b) cell in the null panel
  <species>/spectra/null_<mode>.zball_sorted.npy   the same for the Hamming-1 ball statistic (degenerate pairs)
Usage: python tmap/null_sim.py <species_dir> [mode=nc] [K=5] [seed=1]
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats import ball_z, spacing_z  # noqa: E402

TMAP = os.path.expanduser('~/.cache/tmap_target/release/tmap')
ENC = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'):
    ENC[ch] = i
    ENC[ch + 32] = i


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


def runs(x):
    """(start, end) of maximal non-N runs of an encoded sequence."""
    ok = np.concatenate([[False], x < 4, [False]])
    d = np.diff(ok.astype(np.int8))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def main():
    d = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'nc'
    K = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    rng = np.random.default_rng([20261012, int(sys.argv[4]) if len(sys.argv) > 4 else 1])
    paths = [l.strip() for l in open(os.path.join(d, f'genomes_{mode}.txt')) if l.strip()]
    # fit: counts of (K-mer context, next base) over both strands of every non-N run
    cnt = np.zeros((4 ** K, 4), np.float64)
    seqs = []
    for p in paths:
        recs = [(n, ENC[np.frombuffer(s.encode(), np.uint8)]) for n, s in read_fasta(p)]
        seqs.append(recs)
        for _, x in recs:
            rcx = np.where(x[::-1] > 3, 4, 3 - x[::-1]).astype(np.uint8)       # reverse complement, N kept
            for xx in (x, rcx):
                for a, b in runs(xx):
                    r = xx[a:b].astype(np.int64)
                    if len(r) <= K:
                        continue
                    ctx = np.zeros(len(r) - K, np.int64)
                    for j in range(K):
                        ctx = ctx * 4 + r[j:len(r) - K + j]
                    np.add.at(cnt, (ctx, r[K:]), 1)
    cnt += 0.5
    cum = np.cumsum(cnt / cnt.sum(axis=1, keepdims=True), axis=1)
    start = cnt.sum(axis=1) / cnt.sum()
    start_cum = np.cumsum(start)
    outd = os.path.join(d, f'null_{mode}')
    os.makedirs(outd, exist_ok=True)
    lst = []
    for p, recs in zip(paths, seqs):
        out = os.path.join(outd, os.path.basename(os.path.dirname(p)) + '.fna')
        with open(out, 'w') as fh:
            for name, x in recs:
                y = np.full(len(x), 4, np.uint8)
                for a, b in runs(x):
                    n = b - a
                    u = rng.random(n + 1)
                    ctx = int(np.searchsorted(start_cum, u[0]))
                    first = [(ctx >> (2 * (K - 1 - j))) & 3 for j in range(K)]
                    seg = np.empty(n, np.uint8)
                    m = min(K, n)
                    seg[:m] = first[:m]
                    mask = 4 ** K - 1
                    for i in range(K, n):
                        b2 = int(np.searchsorted(cum[ctx], u[i]))
                        b2 = min(b2, 3)
                        seg[i] = b2
                        ctx = ((ctx << 2) | b2) & mask
                    y[a:b] = seg
                fh.write(f'>{name}\n' + ''.join('ACGTN'[v] for v in y) + '\n')
        lst.append(out)
    open(os.path.join(outd, 'genomes.txt'), 'w').write('\n'.join(lst) + '\n')
    subprocess.run([TMAP, 'spectrum', os.path.join(outd, 'genomes.txt'), os.path.join(outd, 'spec'), '--core', '4',
                    '--gmax', '96'], check=True)
    P = np.fromfile(os.path.join(outd, 'spec.pairs.u32'), np.uint32).reshape(97, 256, 256)
    np.save(os.path.join(d, 'spectra', f'null_{mode}.zsorted.npy'), np.sort(spacing_z(P).ravel()).astype(np.float32))
    np.save(os.path.join(d, 'spectra', f'null_{mode}.zball_sorted.npy'), np.sort(ball_z(P).ravel()).astype(np.float32))
    print(f'null {mode}: K={K}, {len(lst)} synthetic strains -> null_{mode}.zsorted.npy')


if __name__ == '__main__':
    main()
