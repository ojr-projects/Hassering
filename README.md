# Hassering

Hassering is a genome browser built on the genome's own statistical tendencies. For a panel of divergent strains of one
species, it maps every exact sparse pattern that recurs more (or with more regular spacing) than an order-5 Markov twin
genome predicts. No windows are used, and no repeat masking. It then draws the result as a circular genome with the
pattern hierarchy inside the circle.

It is an instrument, not a classifier: it shows exact structure with E-values and applies no pass/fail gates. The
design and the literature behind it are in `docs/GENOME_TENDENCY_MAP_REVIEW.md` (Revisions 1 and 2 are current).

## What it computes

- **Domain:**
  - **non-coding (default):** CDS from Pyrodigal and RNA genes from Bakta are masked to N; coding hits are mapped back
    afterwards;
  - **whole-genome:** shown alongside. Codon periodicity dominates it, which is why non-coding is the default.
- **Spectrum** (`rust/tmap`): exact counts of gapped cores (k = 4–5 within span 13), plus an offset-pair spectrum of
  4-mer cores at gaps 0–96. Both strands, anchored at every base.
- **Statistics** (`tmap/stats.py`):
  - spacing z with a placement-normalised baseline;
  - a Hamming-1 "degenerate ball" statistic;
  - gap-stratified E-values from an order-5 Markov twin panel (`tmap/null_sim.py`).
- **Composites** (`tmap/composites.py`): seeds grow by an exact Apriori chain. The *broad* level holds ≥ 90% of seed
  occurrences; *strict* is the full chain. Families sharing loci merge into *structures*. Degenerate (Hamming-1) groups
  are also formed.
- **Cross-strain bundles** (`tmap/bundles.py`): orthologous occurrences are linked by same-gene adjacency at a similar
  distance, then by identical surrounding pattern structure, then by MinHash of the ±250 bp sequence.
- **Known-motif library** (`tmap/motif_library.json`, built by `tmap/build_library.py`): 486 motifs (PRODORIC 354,
  RegulonDB 95, classic elements 37). Nearest matches are calibrated against column-shuffled queries.

## Viewer

```bash
python tmap/serve.py [port=8765] [panels_root=~/data/tmap/panels]     # then open http://localhost:8765/
```

The viewer has three main views:
- **Circle:** strains as rings, gene arrows, and a backbone that resolves to bases at deep zoom with motif boxes.
  Pointer-anchored wheel zoom.
- **Family list:** filtered by E-value, sortable by E, bits, length or occurrences. Each family has a sequence logo,
  per-strain density bars, and a copy-able structure report.
- **Aligned view:** occurrences one row per locus, oriented to the structure's lead family, grouped in cross-strain
  bundles, with Bakta genes and other motifs shown as blocks.

`tmap/build_viewer.py` writes a standalone HTML file instead (occurrences capped at 4,000 per family).

## Pipeline

```bash
python tmap/panel_fetch.py ~/data/tmap/panels            # NCBI non-MAG assemblies per species (list in the script)
python tmap/panel_diverse.py ~/data/tmap/panels "<Species>" 15   # MinHash farthest-point panel for species with many genomes
tmap/annotate_panels.sh                                   # Bakta (light DB) for every genome; SPECIES_GLOB to restrict
(cd rust/tmap && CARGO_TARGET_DIR=~/.cache/tmap_target cargo build --release)
tmap/run_species.sh <Species_dir> ["nc all"] [top=150]    # mask, spectra, twin null, composites, viewer data
tmap/refresh_views.sh                                     # rebuild every species' view JSON (after viewer_data changes)
```

**Environments** (micromamba):
- `nanomotif`: Python with numpy and the pipeline;
- `bakta`: Bakta ≥ 1.10 with AMRFinderPlus, plus pyrodigal for `mask.py`;
- `jscheck`: node + jsdom, for viewer smoke tests.

## Layout

| path | what |
|---|---|
| `rust/tmap/` | `tmap spectrum <genomes.txt> <prefix> [--core 4] [--gmax 96] [--threads N]`: core counts and offset-pair spectrum |
| `tmap/` | panel fetching, masking, statistics, composites, bundles, view JSON, the server, the viewer template |
| `docs/` | the design review and general development insights |
