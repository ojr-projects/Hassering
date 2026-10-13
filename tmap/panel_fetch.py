# %% panel_fetch.py: multi-strain species panels for the tendency map (docs/GENOME_TENDENCY_MAP_REVIEW.md, Rev. 1)
"""
For each species: list NCBI assemblies (no MAGs, no atypical), rank them (complete genome first, then fewest contigs),
keep up to N with <= MAX_CONTIGS contigs, and download FASTA + GFF3 into <root>/<species>/<accession>/.
A panel.tsv records accession, strain, level, contigs, N50, size and GC. Dereplication by Mash/ANI follows in a later
step; this step only fetches.
Usage: python tmap/panel_fetch.py <root> [N=15] [MAX_CONTIGS=100]
"""
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import zipfile

SPECIES = ['Bifidobacterium adolescentis', 'Barnesiella intestinihominis', 'Bacteroides fragilis', 'Phocaeicola vulgatus',
           'Collinsella aerofaciens', 'Fusicatenibacter saccharivorans', 'Blautia wexlerae', 'Oliverpabstia intestinalis',
           'Lachnospira rogosae']
API = 'https://api.ncbi.nlm.nih.gov/datasets/v2'


def get(url, binary=False, tries=4):
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.read() if binary else json.load(r)
        except Exception:  # noqa: BLE001
            time.sleep(3 * (t + 1))
    raise RuntimeError(url)


def reports(species):
    q = urllib.parse.quote(species)
    out, token = [], None
    while True:
        url = (f'{API}/genome/taxon/{q}/dataset_report?filters.is_metagenome_derived=METAGENOME_DERIVED_EXCLUDE'
               f'&filters.exclude_atypical=true&page_size=1000' + (f'&page_token={token}' if token else ''))
        r = get(url)
        out += r.get('reports', [])
        token = r.get('next_page_token')
        if not token:
            return out


def main():
    root = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    maxc = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    for sp in SPECIES:
        reps = reports(sp)
        rows, seen = [], set()
        for g in reps:
            a, s = g.get('assembly_info', {}), g.get('assembly_stats', {})
            key = g['accession'].split('_')[1].split('.')[0]          # GCA/GCF pairs are one assembly
            if key in seen:
                continue
            seen.add(key)
            rows.append(dict(acc=g['accession'], strain=g.get('organism', {}).get('infraspecific_names', {}).get('strain', ''),
                             level=a.get('assembly_level', ''), contigs=int(s.get('number_of_contigs', 10 ** 6)),
                             n50=int(s.get('contig_n50', 0)), size=int(s.get('total_sequence_length', 0)),
                             gc=s.get('gc_percent', '')))
        rows.sort(key=lambda r: (r['level'] != 'Complete Genome', r['contigs'], -r['n50']))
        keep = [r for r in rows if r['contigs'] <= maxc][:n]
        d = os.path.join(root, sp.replace(' ', '_'))
        os.makedirs(d, exist_ok=True)
        for r in keep:
            gd = os.path.join(d, r['acc'])
            if os.path.exists(os.path.join(gd, 'genome.fna')):
                continue
            z = get(f'{API}/genome/accession/{r["acc"]}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF',
                    binary=True)
            os.makedirs(gd, exist_ok=True)
            with zipfile.ZipFile(io.BytesIO(z)) as zf:
                for m in zf.namelist():
                    if m.endswith('.fna'):
                        open(os.path.join(gd, 'genome.fna'), 'wb').write(zf.read(m))
                    elif m.endswith('.gff'):
                        open(os.path.join(gd, 'genomic.gff'), 'wb').write(zf.read(m))
        with open(os.path.join(d, 'panel.tsv'), 'w') as fh:
            fh.write('acc\tstrain\tlevel\tcontigs\tn50\tsize\tgc\n')
            for r in keep:
                fh.write('\t'.join(str(r[k]) for k in ('acc', 'strain', 'level', 'contigs', 'n50', 'size', 'gc')) + '\n')
        lv = sum(r['level'] == 'Complete Genome' for r in keep)
        print(f'{sp}: {len(rows)} non-MAG assemblies; panel {len(keep)} ({lv} complete, contigs '
              f'{min([r["contigs"] for r in keep], default=0)}-{max([r["contigs"] for r in keep], default=0)})', flush=True)


if __name__ == '__main__':
    main()
