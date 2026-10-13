#!/bin/bash
# Full tendency-map pipeline for one species panel: masks, spectra, Markov twin null, composites, viewer data.
# Usage: tmap/run_species.sh <species_dir_name> [modes="nc all"] [top=150]
set -eu
cd "$(dirname "$0")/.."
SP=$1; MODES=${2:-"nc all"}; TOP=${3:-150}
D=$HOME/data/tmap/panels/$SP
PY=$HOME/micromamba/envs/nanomotif/bin/python
TM=$HOME/.cache/tmap_target/release/tmap
$HOME/micromamba/envs/bakta/bin/python tmap/mask.py $HOME/data/tmap/panels $SP
mkdir -p $D/spectra
for m in $MODES; do
  $TM spectrum $D/genomes_$m.txt $D/spectra/$m --core 4 --gmax 96 2>&1 | tail -1
  $PY tmap/null_sim.py $D $m 5 1
  $PY - <<PYEOF
import sys, numpy as np; sys.path.insert(0, 'tmap'); import stats
P = np.fromfile('$D/null_$m/spec.pairs.u32', np.uint32).reshape(97, 256, 256)
np.save('$D/spectra/null_$m.zgap.npy', stats.null_by_gap(P)); np.save('$D/spectra/null_$m.zballgap.npy', stats.null_by_gap(P, True))
PYEOF
  $PY tmap/composites.py $D $m $TOP
  $PY tmap/viewer_data.py $D $m
done
echo "done $SP"
