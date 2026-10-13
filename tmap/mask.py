# %% mask.py: coding/non-coding domains per panel genome (docs/GENOME_TENDENCY_MAP_REVIEW.md, Revision 2B)
"""
For every genome of a species panel: call protein-coding genes with Pyrodigal (single-genome mode, trained on that
genome; the same CDS caller Bakta uses), then write
  cds.bed         CDS intervals (record, start, end, strand)
  genome_nc.fna   the genome with every CDS base and every RNA-gene base (tRNA, rRNA, ncRNA, tmRNA from Bakta, when
                  annotated) replaced by N: the non-coding, non-gene sequence domain
  mask.done       which mask was applied ('cds+rna' or 'cds'); a changed definition rebuilds genome_nc.fna
plus <species>/genomes_all.txt and genomes_nc.txt, the input lists for `tmap spectrum`.
Bakta's full annotation (RNA genes, CRISPR, oriC, products) is kept separately for the browser layer.
Usage: python tmap/mask.py <panel_root> <species_dir_name>
"""
import os
import sys

import pyrodigal


def read_fasta(p):
    recs, name, buf = [], None, []
    for line in open(p):
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


RNA_TYPES = {'tRNA', 'rRNA', 'ncRNA', 'tmRNA'}


def rna_genes(gd):
    """RNA gene intervals from Bakta's GFF3 (None if Bakta has not annotated this genome)."""
    p = os.path.join(gd, 'bakta', 'genome.gff3')
    if not os.path.exists(p):
        return None
    out = []
    for line in open(p):
        if line.startswith('##FASTA'):
            break
        f = line.split('\t')
        if len(f) > 4 and f[2] in RNA_TYPES:
            out.append((f[0], int(f[3]) - 1, int(f[4])))
    return out


def main():
    root, sp = sys.argv[1], sys.argv[2]
    d = os.path.join(root, sp)
    gl_all, gl_nc = [], []
    for a in sorted(os.listdir(d)):
        gd = os.path.join(d, a)
        if not (a.startswith('GC') and os.path.exists(os.path.join(gd, 'genome.fna'))):
            continue
        recs = read_fasta(os.path.join(gd, 'genome.fna'))
        nc_path = os.path.join(gd, 'genome_nc.fna')
        rna = rna_genes(gd)
        stamp = os.path.join(gd, 'mask.done')
        want = 'cds+rna' if rna is not None else 'cds'
        if os.path.exists(nc_path) and (not os.path.exists(stamp) or open(stamp).read().strip() != want):
            os.remove(nc_path)                         # rebuild when the mask definition changes
        if not os.path.exists(nc_path):
            finder = pyrodigal.GeneFinder()
            finder.train(*[s.encode() for _, s in recs if len(s) >= 20000][:1] or [recs[0][1].encode()])
            with open(os.path.join(gd, 'cds.bed'), 'w') as fb, open(nc_path + '.tmp', 'w') as fo:
                ncod, ntot = 0, 0
                for name, s in recs:
                    m = bytearray(s.encode())
                    for g in finder.find_genes(s.encode()):
                        a0, b0 = g.begin - 1, g.end
                        fb.write(f'{name}\t{a0}\t{b0}\t{"+" if g.strand == 1 else "-"}\n')
                        m[a0:b0] = b'N' * (b0 - a0)
                    for r, a0, b0 in (rna or []):          # RNA genes (Bakta) are masked too: protein and RNA genes out
                        if r == name:
                            m[a0:b0] = b'N' * (b0 - a0)
                    ncod += sum(1 for ch in m if ch == 78)
                    ntot += len(m)
                    fo.write(f'>{name}\n{m.decode()}\n')
            os.replace(nc_path + '.tmp', nc_path)
            open(stamp, 'w').write(want + '\n')
        gl_all.append(os.path.join(gd, 'genome.fna'))
        gl_nc.append(nc_path)
    open(os.path.join(d, 'genomes_all.txt'), 'w').write('\n'.join(gl_all) + '\n')
    open(os.path.join(d, 'genomes_nc.txt'), 'w').write('\n'.join(gl_nc) + '\n')
    print(f'{sp}: {len(gl_all)} genomes masked')


if __name__ == '__main__':
    main()
