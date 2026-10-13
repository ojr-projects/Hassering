# %% stats.py: spacing statistics, Hamming-ball (degenerate) statistics, E-values and IUPAC for the tendency map
"""
spacing_z(P)      z = (count - base) / sqrt(base + 1), base = the cell's median rate per valid placement over gaps x the
                  placements at g (placements shrink with gap in the non-coding domain)
ball_z(P)         the same statistic for the Hamming-1 ball of both cores (a' within 1 mismatch of a, b' of b, same g):
                  degenerate pairs; baseline = sum of each ball cell's own placement-scaled baseline (separable ball operator)
e_value(z, null)  expected number of null cells with z' >= z (null = sorted null z of a Markov-simulated twin panel of the
                  same size); beyond the null maximum, an exponential fit to the top 1,000 null cells extrapolates
iupac_ball(P, g, a, b)  per core position, the bases carrying >= 25% of the excess counts over the ball -> IUPAC cores
"""
import numpy as np

C = 4
IUPAC = {frozenset('A'): 'A', frozenset('C'): 'C', frozenset('G'): 'G', frozenset('T'): 'T', frozenset('AG'): 'R',
         frozenset('CT'): 'Y', frozenset('GC'): 'S', frozenset('AT'): 'W', frozenset('GT'): 'K', frozenset('AC'): 'M',
         frozenset('CGT'): 'B', frozenset('AGT'): 'D', frozenset('ACT'): 'H', frozenset('ACG'): 'V', frozenset('ACGT'): 'N'}


def placements(P):
    """Valid placements per gap (all pairs counted at that gap): longer footprints fit fewer N-free stretches."""
    return P.reshape(P.shape[0], -1).sum(axis=1).astype(np.float64)


def spacing_base(P):
    """Expected count per (g, a, b): the cell's median rate per placement over gaps, times the placements at g."""
    Pf = P.astype(np.float64)
    N = placements(Pf)
    rate = Pf / N[:, None, None]
    return np.median(rate, axis=0)[None] * N[:, None, None]


def spacing_z(P):
    base = spacing_base(P)
    return (P.astype(np.float64) - base) / np.sqrt(base + 1.0)


def ball_axis(X, axis):
    """Hamming-1 ball sum of a 4-mer axis (self included): sum over j of the 4 variants at position j, minus 3 x self."""
    X = np.moveaxis(X, axis, -1)
    sh = X.shape
    V = X.reshape(sh[:-1] + (4,) * C)
    out = -(C - 1) * V
    for j in range(C):
        out = out + V.sum(axis=len(sh) - 1 + j, keepdims=True)
    return np.moveaxis(out.reshape(sh), -1, axis)


def ball_counts(P):
    return ball_axis(ball_axis(P.astype(np.float64), 1), 2)


def ball_z(P):
    """Degenerate statistic: (sum of ball counts - sum of each ball cell's own median over gaps) / sqrt(that + 1)."""
    Pf = P.astype(np.float64)
    B = ball_counts(Pf)
    base = ball_counts(spacing_base(Pf))
    return (B - base) / np.sqrt(base + 1.0)


def e_value(z, null_sorted, n_obs_cells=None):
    """Expected number of cells with statistic >= z in a null panel of the same size (two orientations of each canonical
    pair are both cells in both panels, so the scale is consistent)."""
    n = len(null_sorted)
    k = n - np.searchsorted(null_sorted, z, side='left')
    if k > 0:
        return float(k)
    top = null_sorted[-1000:].astype(np.float64)
    ranks = np.arange(1000, 0, -1, dtype=np.float64)              # number of null cells >= each top value
    slope, icpt = np.polyfit(top, np.log(ranks), 1)
    return float(np.exp(icpt + slope * z)) if slope < 0 else 0.0


def dec(x):
    return ''.join('ACGT'[(x >> (2 * (C - 1 - i))) & 3] for i in range(C))


def neighbours(x):
    out = [x]
    for j in range(C):
        sh = 2 * (C - 1 - j)
        cur = (x >> sh) & 3
        for b in range(4):
            if b != cur:
                out.append((x & ~(3 << sh)) | (b << sh))
    return out


def iupac_ball(P, base, g, a, b, share=0.25):     # base: spacing_base(P)
    """IUPAC forms of cores a and b from the excess counts of their Hamming-1 variants at gap g."""
    wa, wb = np.zeros((C, 4)), np.zeros((C, 4))
    for a2 in neighbours(a):
        ex = max(P[g, a2, b] - base[a2, b], 0.0)            # base: the gap-g slice of spacing_base
        for j in range(C):
            wa[j, (a2 >> (2 * (C - 1 - j))) & 3] += ex
    for b2 in neighbours(b):
        ex = max(P[g, a, b2] - base[a, b2], 0.0)
        for j in range(C):
            wb[j, (b2 >> (2 * (C - 1 - j))) & 3] += ex

    def form(w, x):
        s = ''
        for j in range(C):
            tot = w[j].sum()
            keep = {'ACGT'[k] for k in range(4) if tot > 0 and w[j, k] / tot >= share}
            keep.add('ACGT'[(x >> (2 * (C - 1 - j))) & 3])
            s += IUPAC[frozenset(keep)]
        return s
    return form(wa, a), form(wb, b)


def null_by_gap(P_null, ball=False):
    """Sorted null statistic per gap: array (G, cells) for gap-stratified p-values."""
    Zn = ball_z(P_null) if ball else spacing_z(P_null)
    return np.sort(Zn.reshape(Zn.shape[0], -1), axis=1).astype(np.float32)


def e_value_gap(z, g, nbg):
    """Gap-stratified E-value: p = share of null cells at the same gap with z' >= z (exponential tail fit to that gap's
    top 500 cells beyond its maximum), times the number of canonical cells tested (half of all cells)."""
    row = nbg[g]
    n = len(row)
    k = n - np.searchsorted(row, z, side='left')
    if k > 0:
        p = k / n
    else:
        top = row[-500:].astype(np.float64)
        ranks = np.arange(500, 0, -1, dtype=np.float64)
        slope, icpt = np.polyfit(top, np.log(ranks), 1)
        p = float(np.exp(icpt + slope * z)) / n if slope < 0 else 0.0
    return p * nbg.size / 2
