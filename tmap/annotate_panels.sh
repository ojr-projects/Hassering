#!/bin/bash
# Uniform Bakta annotation of every panel genome (docs/GENOME_TENDENCY_MAP_REVIEW.md, Revision 2B).
# Output: <genome_dir>/bakta/genome.{gff3,tsv,...}. Waits for the light DB; 6 genomes in parallel, 2 threads each.
set -u
DB=$HOME/data/tmap/bakta_db/db-light
B=$HOME/micromamba/envs/bakta/bin
until [ -f $DB/version.json ]; do sleep 30; done
export PATH=$B:$PATH
ls -d $HOME/data/tmap/panels/${SPECIES_GLOB:-*}/GC*/ | xargs -P 6 -I{} sh -c '
  d={}; [ -f $d/bakta/genome.gff3 ] && exit 0
  bakta --db '"$DB"' --output $d/bakta --prefix genome --threads 2 --skip-plot --keep-contig-headers --force $d/genome.fna \
    > $d/bakta.log 2>&1 && echo "done $d" || echo "FAILED $d"'
touch $HOME/data/tmap/panels/annotate_${SPECIES_TAG:-all}.done
