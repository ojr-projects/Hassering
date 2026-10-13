#!/bin/bash
# Wait for running species pipelines to finish, then rebuild every species' viewer data (both modes).
while pgrep -f "[t]map/run_species.sh" >/dev/null; do sleep 60; done
for d in $HOME/data/tmap/panels/*/; do for m in nc all; do
  [ -f $d/spectra/$m.families.json ] && $HOME/micromamba/envs/nanomotif/bin/python tmap/viewer_data.py $d $m
done; done
echo refreshed
