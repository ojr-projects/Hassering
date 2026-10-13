# Genome tendency map: a self-referential, hierarchical map of sparse-pattern tendencies

*Literature and methods review plus a design proposal, 2026-10-12.*

**The idea, in one paragraph (user brief).** An all-vs-all comparison of a genome against itself, Comsequent-style
rather than BLAST-style. Find every region that is *cross-reactive* with others because they share loose,
statistical, sparse-pattern tendencies spread over tens of bases, not because they are homologous. Break the genome
down progressively into groups defined by those tendencies, and build a hierarchical map of region types, de novo.
Do not mask repeats; they should fall out as their own groups early. Start with microbes, and show the result as a
circular genome with the discovered hierarchy drawn inside the circle: a genome browser that reports what the
sequence itself "says" before any annotation.

**Verification legend.** Claims are from sources fetched during this review unless marked:
- **[bg]** from background knowledge; verify before relying on it.
- **[derived]** worked out here.

The literature search was a quick pass (about 15 queries), not a systematic review. Absence of prior art below means
"not found in that pass".



> **Revision 2 (2026-10-12, user direction). It supersedes §4.1–§4.3 (windows) and adds the non-coding mode.**
>
> **A. No windows. Structure is anchored at every base and defined by exact geometry.**
> - **Primitive:** an exact occurrence of a sparse pattern at (strain, contig, position, strand). Every base anchors
>   every pattern that starts there. There are no units and no boundaries.
> - **Level 0 (cores):** all gapped patterns (k ≤ 5 in span ≤ 13; contiguous to L = 8), with exact counts and exact
>   surprise against the strain-pooled context model of the active sequence domain.
> - **Level 1 (the tens-of-bases scale):** the exact co-occurrence spectrum of pairs of cores at a fixed relative
>   offset d ≤ D (≈ 100 bp) and a fixed orientation. Example: "TTGA, then 17 bp later, TATA".
>   - This reaches long, sparse spans without enumerating C(60, 6)·4⁶ patterns.
>   - The expectation comes from the two cores' own counts and the local context (the independence/IPF residual).
>   - Counting uses sorted occurrence lists of the selected cores, so it is exact. Densely, 3-mer cores give 64² ×
>     2D cells.
> - **Level n:** composites extend by adding one core at a fixed offset (Apriori-style; a child's occurrences are an
>   exact subset of its parent's).
> - **Hierarchy:** pattern inclusion, giving an exact lattice. Cross-reactivity edges join composites that occupy
>   the same loci, and the same orthologous loci across strains.
> - **Regions emerge:** footprints of composite occurrences. Overlapping footprints of related nodes merge into
>   tracts whose boundaries come from the data.
> - **D is a reach limit, not a boundary.**
>
> **B. Sequence domain: non-coding mode (a default candidate), with coding mapped back.**
> - **Annotation:** every panel genome is annotated uniformly with Bakta 1.12 (light DB v6). NCBI GFFs are missing
>   for 56 of 116 genomes, and mixed annotation sources would make coding boundaries inconsistent across strains.
>   Bakta also provides tRNA/rRNA/ncRNA, CRISPR arrays and oriC for the browser layer.
> - **Non-coding domain:** the complement of protein-coding CDS. RNA genes are included and flagged.
>   - A pattern or composite counts only if its whole footprint lies inside the domain.
>   - The context model is fitted on the domain only.
> - **Mapping back:** each non-coding node's patterns are counted exactly in CDS too, reporting density in coding vs
>   non-coding, strand and codon phase. "This tendency is intergenic-specific" vs "it is everywhere" is then visible,
>   with no corruption of the discovery step by codon structure.
> - **Both modes are computed at first** (non-coding, and whole genome) so the two maps can be compared.

> **Revision 1 (2026-10-12, user direction). It supersedes §7–§8 and sharpens §4. Read it first.**
>
> 1. **This is an instrument, not a classifier.** There are no pass/fail gates. The goal is to show **raw, exact
>    structure**: real self-interactive tendencies in the sequence, often with unknown explanations. Those tendencies
>    are not necessarily distinct "element types", and their separation from other structure is variable. Statistics
>    say how surprising a structure is; they do not decide whether it is shown. Calibration (§4.2 step 5) becomes an
>    annotation on each node ("survives a context-matched null at level x"), not a stopping gate.
> 2. **Exact first; fuzzy steps are optional.** The default node definition is exact (§4.7):
>    - groups of units defined by **exact shared pattern sets** (closed sets of the unit × pattern occurrence
>      relation);
>    - ranked by exact surprise against the context model.
>    Mixture, NMF and topic splits (§4.2 step 4, §4.4) remain optional alternative views.
> 3. **A multi-strain species is the native default unit (§4.8).** Several divergent complete genomes of one species
>    are analysed together:
>    - the context model and the pattern statistics are pooled across strains, so one individual's variation is not
>      overfitted;
>    - every node reports how consistently it occurs across strains: core, accessory or strain-private.
> 4. **Organisms (replaces §7 and §8).** Gut microbes first, not *E. coli*. Complete-genome availability
>    (NCBI Datasets, 2026-10-12; "all" includes drafts and MAGs):
>
>    | species | complete genomes | all assemblies | why |
>    |---|---|---|---|
>    | *Bifidobacterium adolescentis* | 157 | 2,199 | user's first choice; abundant; plenty of complete strains |
>    | *Barnesiella intestinihominis* | 6 | 844 | user's choice; few complete, so strain panel ≈ 6 + best drafts |
>    | *Bacteroides fragilis* | 232 | 3,499 | **self-recombination**: many invertible promoters ("invertons"; Jiang et al. 2019, *Science* 363:181; PhaseFinder found 4,686 invertons in 54,875 genomes, enriched in host-associated bacteria). Inverted repeat pairs flanking promoters are exactly the self-structure the map should show |
>    | *Phocaeicola vulgatus* | 98 | 3,336 | Bacteroidales with invertons; second self-recombination panel |
>    | *Collinsella aerofaciens* ("the monster") | 18 | 758 | strain heterogeneity near the species boundary (ANI ~95%, DDH < 70% for some strains; [MDPI 2019](https://www.doi.org/10.3390/MICROORGANISMS7030078)); later |
>    | *Faecalibacterium prausnitzii* | 51 | 2,033 | abundant; known species-complex structure (optional) |
>
>    **Aristaeellaceae** (GTDB family, order Christensenellales; checked 2026-10-12): **dropped for now.**
>    - The family has 1,705 genomes in 590 species clusters.
>    - The abundant human-gut cluster, UBA11524 sp000437595 ("*Faecalibacterium* sp. CAG:74"), has 265 genomes,
>      all MAGs: 10 with ≤ 20 contigs and none closed.
>    - The only closed isolates are rumen *Aristaeella* (*A. lactis*, *A. hokkaidonensis*, *A.* sp900322155), with at
>      most 2 per species.
>    - Revisit when closed long-read MAGs or isolates of UBA11524 appear.
>
>    **Also poorly understood:** bifidobacterial prophages carry a putative tyrosine-recombinase shufflon (Rin) that
>    diversifies receptor-binding proteins ([Mavrich et al. 2018](https://cora.ucc.ie/items/7ffe42b8-7ef5-4aed-8934-0b026e270410)).
>    That is a self-recombination structure inside the *Bifidobacterium* panel itself.

---

## 1. Summary

1. **This is buildable now, at microbial scale, on the workstation.**
   - The data structure exists: exact gapped-pattern counting at ~2×10⁹ increments per core-second (HASTUR's Rust
     core).
   - The statistics exist: HASTUR's permutation max-T calibration.
   - A 2–10 Mb genome needs seconds to minutes per hierarchy level [derived, §5].
2. **The reframing that makes it work where Comsequent failed.**
   - Comsequent compared a set against a background, and failed because background structure saturated its tests:
     every pattern became significant under any composition mismatch (`docs/final_reflections.md` §9).
   - A self-map tests **each node against its parent**, never against the whole genome. Dominant structure (coding
     vs non-coding, codon periodicity, GC/isochore, strand skew, repeat families) is explained at the top levels by
     construction. Deeper levels are conditional residuals.
   - The background problem becomes the top of the map instead of a confound.
3. **HASTUR's one proven asset, calibration, is exactly what a self-map lacks elsewhere.**
   - Every clustering method will split forever.
   - A permutation-calibrated split test turns "how deep should the tree go" into a family-wise-controlled decision.
   - The output is a map whose every node is a statistically supported region type, defined by an explicit,
     readable dictionary of sparse patterns.
4. **Prior art is adjacent, not identical.**
   - Composition-based fragment clustering (BLSOM, metagenome binning) works at kilobase scale on 4–5-mer frequencies.
   - Entropic recursive segmentation finds contiguous compositional domains.
   - Genomic-LM embeddings and sparse-autoencoder features learn region types opaquely, and the deep-model review
     shows they are mostly k-mer statistics.
   - Not found: a **calibrated, hierarchical, parent-conditional, gapped-pattern map at tens-of-bases resolution with
     an interpretable dictionary per node**.
5. **The main risk is not compute but meaning.** A map that rediscovers coding/non-coding, rRNA, REP elements and
   promoters is a correct but unsurprising browser. The value rests on whether the deeper nodes are (a) real, which
   calibration checks, and (b) informative, which annotation overlap and cross-species recurrence check.
   **§7 states the cheapest decisive test.**

---

## 2. Prior art

### 2.1 Composition signatures of genome fragments (kilobase scale)

- **BLSOM: batch-learning self-organising maps of oligonucleotide frequencies.**
  - Clusters 1–100 kb fragments by 4–5-mer frequencies (PCA-initialised, order-independent), recovering phylotype
    without labels.
  - Scale: ~10⁶ fragments from 101 eukaryotes; >5×10⁶ microbial fragments
    ([PMC3967822](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3967822/);
    [DSJ 2021](https://datascience.codata.org/article/10.5334/dsj-2021-029/)).
  - **Relevance:** the "genome signature" idea at coarse scale. It sees species and islands, not elements tens of
    bases long.
- **Metagenome binning.**
  - VAMB is a VAE over tetranucleotide frequencies plus coverage, then clustering
    ([VAMB](https://communities.springernature.com/posts/improved-metagenome-binning-and-assembly-using-deep-variational-autoencoders)).
    AAMB is an adversarial variant
    ([bioRxiv 2023](https://www.biorxiv.org/content/10.1101/2023.02.27.527078.full.pdf)).
  - The VAMB authors note that VAE assumptions are "not very amenable to clustering in general".
  - **Relevance:** proven engineering for clustering short-sequence composition. They operate on contigs (kb), and
    composition is the signal rather than the nuisance.

### 2.2 Segmentation of a genome into compositional domains

- **Jensen–Shannon recursive segmentation** (Azad, Ramaswamy et al.).
  - Splits a sequence recursively into compositionally homogeneous domains.
  - Higher-order Markov versions find genomic island and prophage boundaries in *E. coli* K-12
    ([Markov models of genome segmentation](https://repository.ias.ac.in/98848)).
  - Domains can "share compositional similarity with one or more distant domains" ([arXiv physics/0207113](https://arxiv.org/pdf/physics/0207113)). That
    is the cross-reactivity idea at domain scale.
  - Arvey et al. 2009 built a top-down island detector on the same principle
    ([PMC2760805](https://pmc.ncbi.nlm.nih.gov/articles/PMC2760805)).
- **Alien_Hunter / IVOM:** interpolated variable-order motif scores across sliding windows, plus a change-point HMM,
  for horizontally transferred regions ([man page](https://manpages.debian.org/jessie/alien-hunter/alien_hunter.1.en.html)).
- **Relevance:**
  - These methods are the right *top level* of the proposed map: islands, prophages and isochores.
  - They are contiguous-segment methods. They do not group scattered short elements by shared tendency, and they
    do not build a hierarchy of types.

### 2.3 Self-similarity without homology search

- **Compression with expert models** (XM; GeCo3).
  - XM blends a background model, a local context model and forward/reverse-complement copy experts. Its authors
    propose it for repeat discovery
    ([XM](https://allisons.org/ll/links/DNAcompress-XM.html); [GeCo3](https://pmc.ncbi.nlm.nih.gov/articles/PMC7657843)).
  - **Relevance:** "region A helps compress region B" is a principled cross-reactivity measure, and it includes
    approximate repeats. It is a candidate *similarity* for the all-vs-all graph (§4.3), and a candidate MDL stopping
    rule.
- **Spaced-word alignment-free comparison** (Morgenstern lab).
  - Uses match/don't-care patterns for distances between whole sequences
    ([spaced words](https://publications.goettingen-research-online.de/ID-32985)).
  - **Relevance:** a theoretical precedent for gapped patterns as a similarity basis. It is used between genomes,
    not within one.
- **Gapped k-mer kernels** (gkm-SVM; LS-GKM; FastSK Monte Carlo approximation).
  - The gkm kernel is an inner product over all positional match subsets in an l-window
    ([Ghandi 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4102394/);
    [FastSK](https://www.biorxiv.org/content/10.1101/2020.04.21.053975v2.full.pdf)).
  - **Relevance:** the natural all-vs-all similarity between windows is exactly this kernel. It is computable through
    counts, since HASTUR already enumerates every gapped cell.
- **De novo repeat finders** (RECON, RED, RepeatScout [bg]; iterative random-PWM dispersed repeats).
  - These find dispersed repeats up to x ≈ 1.5 substitutions per site. REP/BIME palindromes (~40 bp) are the largest
    repeat class in *E. coli* ([PMC3333891](https://pmc.ncbi.nlm.nih.gov/articles/PMC3333891)).
  - **Relevance:** homologous families. In the proposed map they should be the first, most overdispersed splits,
    which makes them a validation layer and not a competitor.

### 2.4 Topic models and mixed membership on k-mers

- LDA on positional k-mers recovers splice-site motifs and intron subtypes
  ([Wang & Mount 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11310860/)).
- LDA on k-mer documents classifies 16S sequences
  ([La Rosa 2015](https://link.springer.com/article/10.1186/1471-2105-16-S6-S2)).
- **Relevance:**
  - A window of 100–200 bp in a bacterium often holds several elements at once: promoter, RBS, start codon, the end
    of an upstream gene's terminator.
  - Mixed membership (a window is a mixture of "tendency topics") is the natural model.
  - A **hierarchical topic model over gapped patterns** is the mixed-membership version of the proposed map (§4.4).

### 2.5 Large models: learned region types

- **GPN** (DNA LM, Brassicales).
  - UMAP of embeddings averaged over 100-bp windows, coloured by gene region; a probe classifies region types
    ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.08.22.504706.full.pdf)).
- **Evo 2 sparse autoencoders.**
  - BatchTopK SAEs on layer 26 give features for prophage regions, CDS, frameshifts, introns and exon boundaries
    ([Goodfire](https://goodfire.ai/blog/interpreting-evo-2); [Evo 2](https://arcinstitute.org/manuscripts/Evo2.pdf)).
  - SAE features beat raw embeddings on 42 of 55 genomic concepts
    ([ICML 2026](https://icml.cc/virtual/2026/poster/62924)).
  - The prokaryotic evidence is mainly the prophage finding. No systematic benchmark on bacterial elements was found
    ([biosecurity probes](https://arxiv.org/pdf/2607.14070)).
- **Bacterial gLMs.**
  - Evo 2's unsupervised representation aligns with supervised CDS calls
    ([GeneLM, Brief Bioinform 2025](https://academic.oup.com/bib/article/26/4/bbaf311/8182749)).
  - The M5 whole-genome bacterial encoder works at single-nucleotide resolution
    ([arXiv 2407.03392](https://arxiv.org/pdf/2407.03392)).
- **What the deep-model review says** (`docs/DNA_deep_models_motif_lit_review.md` §6, §7.5):
  - gLM representations are dominated by k-mer composition, GC and repeats.
  - Masked LMs see a motif only through the mutual information inside it, so short low-information motifs are nearly
    invisible.
  - Transformers approximate variable-order Markov (CTW) estimators.
- **Implication:**
  - SAE features on Evo 2 are the closest existing thing to "an internal representative map of tendencies".
  - It is learned, GPU-heavy and not calibrated, and it offers no statement that a feature is more than chance
    structure.
  - The proposed map computes the same kind of statistics (variable-order context, pattern co-occurrence) explicitly,
    exactly and with calibrated splits.
  - Evo 2 SAE features, or simpler averaged gLM embeddings, are the right **comparator** (§7), not a component.

### 2.6 Known genome-wide "tendencies" in bacteria (validation targets)

These are short, scattered, statistically over-represented signals that a self-map *must* recover:

- **Chi sites:** about 8 nt, the most over-represented oligonucleotide in *E. coli*, with strand bias
  ([Chi distribution](https://researchportal.ip-paris.fr/en/publications/characteristics-of-chi-distribution-on-different-bacterial-genome/)).
- **DNA uptake sequences:** *H. influenzae* has about 1,471 copies of AAGTGCGGT
  ([PMC384778](https://pmc.ncbi.nlm.nih.gov/articles/PMC384778)).
- **REP/BIME palindromes:** about 500 BIMEs and more than 4,000 REPs in MG1655
  ([PMC3333891](https://pmc.ncbi.nlm.nih.gov/articles/PMC3333891)).
- **Other classes:** σ70 promoters (−35/−10 with spacer), Shine–Dalgarno, rho-independent terminators, KOPS/matS/DnaA
  boxes [bg].
- **Cyanobacteria:** the HIP1 octamer (GCGATCGC) is hugely over-represented [bg]. *Nostoc* PCC 7120 and
  *Synechococcus* PCC 7942 are already among the local genomes.
- **Giant viruses:** the Mimivirus early promoter motif AAAATTGA [bg] is a ready truth for a "novel organism" test.

---

## 3. Methods review: what each prior line contributes

| Line | Proven asset | Proven failure | What carries into the map |
|---|---|---|---|
| **Comsequent** (`docs/COMSEQUENT_METHOD.md`) | Exact enumeration of every sparse pattern (k ≤ 7, span ≤ 15); fast (*E. coli* 5/13 in 4 s); lattice algebra; θ (log-linear interaction) sees whole motifs where Möbius I cannot | Contrast tests saturate under composition mismatch; Möbius I is pointwise co-information, not "excess beyond sub-patterns"; exact patterns are diluted on PWM-like signals | The counting engine and the lattice. **θ, or max-ent (IPF) residuals, not Möbius I**, for "tendency beyond sub-patterns". Within-genome conditioning (the methylation case was its one exchangeable win) |
| **HASTUR** (F2, Rust) | Permutation max-T calibration that holds on every null tried; PEXT-coded gapped counting at ~10⁹/s/core; subtraction trick (count one side, derive the other); EM/frame PWM representation | Sensitivity at small N; composition-shifted nulls need matched backgrounds (Mode A) | **The split test**: is this node's best split better than any split of a null that keeps the node's own context statistics? The subtraction trick: count the child, derive the sibling. Representation of each node's defining patterns as logos |
| **PROMETHEUS** | Genome covariate index (GC, CpG o/e, repeat share, rmsk class, TSS distance); matched sets; stratified permutation | Small real-data gain | Covariate tracks as **explanatory overlays** on the map, e.g. "this node is explained by GC" vs "survives GC strata" |
| **Large models** | Learn region types and variable-order context; SAE features are interpretable-ish | Mostly k-mer statistics; uncalibrated; GPU cost; opaque | The *comparator*, plus the insight that variable-order context (CTW/BaMM back-off) is the right null at each node |

---

## 4. Proposed method: the calibrated tendency tree

### 4.1 Units and features
- **Units:** stranded windows (window, strand) of width W, about 100–200 bp, tiled at step W for statistics. Overlap
  is used only for display.
  - **Why stranded:** promoters, RBS, terminators, Chi and replication-oriented motifs are directional. The tree
    should discover orientation, not have it pooled away.
  - Both strands of a window are distinct units. Reverse-complement pairing is recorded so it can be reported.
- **Features:** HASTUR's gapped lattice (k = 4–5 fixed bases in span 13; ~5.6×10⁵ cells), plus contiguous L-mers up
  to 8, counted per unit through the same PEXT window codes. The vocabulary is never materialised per unit; see §5.

### 4.2 Divisive, parent-conditional splitting (one node)
For node N, a set of units:
1. **Context model of N.** Fit a variable-order Markov model (BaMM/CTW back-off, order ≤ 5) to N's pooled sequence.
   It is the null for everything inside N, and it absorbs N's composition.
2. **Tendency statistic per pattern.** The **over-dispersion** of each pattern's per-unit counts across N's units,
   relative to the count variance expected under N's context model.
   - A pattern present at a uniform rate is not a tendency, even if it is enriched genome-wide.
   - A pattern concentrated in a subset of units is a tendency.
   - It is computed as Σ_u c_u(p)² against its expectation, streamed (§5).
   - Optionally use the **IPF residual** (excess beyond all sub-patterns; deep-model review §7.7) so that a 5-base
     tendency is not just its 4-base parent echoed.
3. **Vocabulary.** The top D (≈ 1–5k) over-dispersed patterns, made non-redundant: a sub-pattern is dropped if its
   dispersion is explained by a selected super-pattern.
4. **Split.** Build a units × D sparse count matrix and run a multinomial mixture (2–K components, EM), or NMF, or
   spectral bisection on the gkm-kernel graph restricted to the vocabulary.
5. **Calibration (HASTUR).**
   - Statistic: the best split's held-out log-likelihood gain, maximised over the vocabulary search and the split.
   - Null: the same whole procedure on R pseudo-nodes. Each has the same unit count and lengths, simulated from N's
     context model, or made by k-mer-preserving shuffles of N's units.
   - Fit a Gumbel to the null maxima; split only at family-wise α (e.g. 0.01 per node, with α spent down the tree).
   - The tree's depth is therefore a calibrated quantity.
6. **Recurse** on the children. A child's context model is refitted, so its tendencies are always relative to its
   parent.

**Why repeats need no masking.** Copies of an IS element, rRNA operons and REP/BIME families give extreme
over-dispersion of long patterns. They are split off within the first few levels as tight nodes, and are labelled as
"homologous families" by a copy-identity check (pairwise identity of member units) [derived].

### 4.3 The all-vs-all cross-reactivity graph (BLAST analogue)
- **Similarity:** for units u and v within a node, the gkm-type kernel over the node's vocabulary, weighted by each
  pattern's inverse expected rate under the context model:
  K(u, v) = Σ_p w_p · (c_u(p) − e_u(p)) · (c_v(p) − e_v(p)).
  Unlike BLAST's "these two are homologous", this says "these two over-express the same sparse patterns relative to
  their shared context".
- **Graph:** a k-nearest-neighbour graph built by sparse matrix products (microbes) or locality-sensitive hashing
  (large genomes). Leiden communities at several resolutions give an *agglomerative* hierarchy, which cross-checks
  the divisive tree.
- **Edges for display:** the top cross-reactive pairs, drawn as bundled edges (§6).
- **Optional alternative similarity:** conditional compression gain, i.e. bits saved coding u with v's context
  expert (XM-style). This catches approximate repeats and long-range structure the pattern kernel misses.

### 4.4 Mixed membership (v1, after the hard tree works)
Windows hold several elements. The mixed-membership form is a **hierarchical topic model over gapped patterns**, in
the style of nested-CRP/hLDA [bg], in which each node is a pattern distribution and each unit a path mixture. It is
heavier. Start with the hard tree, then add soft assignment where units sit near a split boundary.

### 4.5 Spatial coupling
Element types have lengths: CDS about 1 kb, promoter regions about 100 bp, REP about 40 bp. A two-pass refinement
re-labels each unit with an HMM along the genome whose states are the tree's leaves. That sharpens boundaries, much
like segmentation literature's change-point HMM (Alien_Hunter) [derived].

### 4.6 Node reports (what the browser shows per node)
- **Defining patterns:** top patterns by enrichment against the parent, with HASTUR EM logos.
- **Positional profile:** position within the window, plus distance to the nearest CDS start/end where annotation
  exists.
- **Orientation:** strand, and position relative to the replication origin/terminus.
- **Membership:** count, GC, copy-identity (homologous family or not), and conservation in related genomes (§8).
- **Annotation overlap:** for validation organisms, overlap with CDS, rRNA/tRNA, REP, promoters, terminators and
  islands.


### 4.7 Exact node definition (default view; Revision 1)
- **The relation:** R(u, p) = 1 if unit u contains pattern p, over the over-dispersed vocabulary of the node, both
  strands, exact matches.
- **Exact groups:** **closed sets**, i.e. formal concepts (U, P). U is exactly the set of units containing all of P,
  and P is exactly the set of vocabulary patterns shared by all of U.
- **The hierarchy:** concepts are ordered by inclusion. The resulting **concept lattice** is the exact,
  heuristic-free structure: no clustering choice, no seed, no tuning.
- **Making it tractable:**
  - closed-itemset mining with a minimum support (e.g. ≥ 5 units, or ≥ 2 strains);
  - a minimum surprise: exact binomial/Poisson for the joint occurrence of P in U against the context model's
    expectation;
  - a maximum pattern-set size.
  - Microbial vocabularies (10³–10⁴ patterns) over 10⁴–10⁵ units are routine for closed-itemset miners
    (LCM-class algorithms [bg]).
- **The tree for display:** the lattice is a DAG. The circle view uses a spanning tree: each concept attaches to its
  most surprising parent. Cross-links are kept as edges.
- **What it shows:** "these N regions (positions, strands) all contain exactly these k sparse patterns, which is
  10^−x surprising under the strain-pooled context". This is the cleanest statement the eye can follow.

### 4.8 Multi-strain species as the default unit (Revision 1)
- **Panel:** S complete genomes (S ≈ 6–20), dereplicated and spread across the species tree. Choose them by Mash
  distance: a farthest-point selection across the ANI range, for example.
- **Pooled context:** the context model (variable-order Markov) is fitted on all strains together. Over-dispersion is
  computed within each strain and combined, so a tendency must be present across strains, not driven by one genome's
  idiosyncrasy.
- **Units:** windows from every strain. Each node records its strain coverage: core (all strains), accessory or
  private.
- **Orthology:** syntenic positions across strains are linked through whole-genome alignment of the panel
  (minimap2/nucmer [bg]). That is used only to *display* "the same locus in other strains", never to define nodes.
- **Display:**
  - one reference strain's circle, with strain-consistency shading per node;
  - optionally a small-multiples ring of circles, one per strain;
  - inversions between strains (the inverton case) appear as positional flips of a node's members, with the
    flanking inverted repeats as their own exact pattern-set node.

---

## 5. Efficiency and scale

The figures below are [derived], from the measured HASTUR counting rate (~2.3×10⁹ gapped-cell increments per core-s)
and E. coli counts (both strands, 9.2×10⁶ positions).

| Step | Cost per tree level | *E. coli* (4.6 Mb), 12-core workstation |
|---|---|---|
| Pooled node counts (all cells) | positions × 715 shapes; subtraction trick for siblings | 6.6×10⁹ → **~0.3 s** |
| Over-dispersion Σ_u c_u(p)² | per unit per shape, a small histogram of ≤ W codes; squared and accumulated | same order: **~1 s** |
| Context-model fit (order 5) | positions | negligible |
| Units × vocabulary matrix (D = 5k) | positions × shapes in vocabulary (~200) | ~2×10⁹ → **~0.2 s** |
| Mixture/EM split | units × nnz × iterations | ~10⁵ units: **seconds** |
| Calibration (R = 20 null nodes) | 20 × the above | **~1 min per level** (most nodes are small) |
| kNN graph (all-vs-all) | sparse products, 9×10⁴ units | **seconds to a minute** |

A tree of about 10 levels and ~100–500 nodes would take **minutes** for a bacterium.

| Genome scale | Units (W = 100, stranded) | Feasibility |
|---|---|---|
| Microbes, 1–10 Mb | 2×10⁴ – 2×10⁵ | **Workstation, minutes.** The intended start. |
| Giant viruses, 1–2.5 Mb | ~5×10⁴ | Workstation; a good novel-organism target |
| Yeasts, 12 Mb | ~2.4×10⁵ | Workstation, ~10 min |
| *Plasmodium*, 23 Mb (13% GC) | ~4.6×10⁵ | Workstation; composition dominates the top levels |
| Fly / worm, 100–140 Mb | ~3×10⁶ | Box; the kNN graph needs locality-sensitive hashing |
| Human, 3 Gb | ~6×10⁷ | Box-scale counting is fine (index precedent); the graph and EM need sketching. Later, if ever |
| Metagenome assemblies | variable | Natural extension: the same tree across many genomes, where nodes become cross-species tendency classes |

**Memory.** Node counts take 5.6×10⁵ cells × 4 B ≈ 2 MB. The units × vocabulary matrix is sparse, about 10⁵ × 10³
nonzeros. Nothing approaches the 15 GB workstation limit at microbial scale.

---

## 6. The visualiser (circular genome, hierarchy inside)

**Layout.**
- The genome is the outer circle, with angle = position, and the origin and terminus marked.
- Annotation rings sit just inside it: CDS by strand, rRNA/tRNA, and GC skew. These are optional and off for a novel
  organism.
- **Hierarchy rings inside.** This is a radial icicle (sunburst-like) plot turned outside-in: ring r is tree depth r,
  with the root innermost or outermost by preference. At each angle, a ring is coloured by the node that the unit at
  that position belongs to at that depth.
  - Coarse levels show long arcs: coding blocks, islands, prophages.
  - Fine levels show speckles: scattered elements.
- **Hierarchical edge bundling** (Holten 2006, *IEEE TVCG* 12(5)).
  - The tree is drawn in the disc, with leaves on the circle at their genome positions.
  - Edges between strongly cross-reactive units are B-splines routed through the tree.
  - REP/BIME families, IS elements and promoter classes then appear as bundles crossing the genome
    ([Holten](https://research.tue.nl/nl/publications/hierarchical-edge-bundles-visualization-of-adjacency-relations-in/)).

**Interaction.**
- Clicking a node highlights its members around the circle.
- A side panel shows defining-pattern logos, the positional profile, counts, orientation and annotation overlap.
- Zoom drops to a linear track at base-pair level, as CGView.js does.

**Implementation.**
- A single-page HTML artifact (D3 v7 with canvas, since 10⁴–10⁵ units are too many for SVG).
- Data comes from a JSON emitted by the tree builder: the node table, per-unit leaf paths, and the top edges.
- Proksee and CGView.js are the circular-map precedent
  ([CGView.js, JOSS](https://www.theoj.org/joss-papers/joss.09930/10.21105.joss.09930.pdf)), but they have no
  in-circle hierarchy. That part would be new.

---

## 7. (Superseded by Revision 1) Original evaluation plan

**Organisms**, all already local:
- *E. coli* K-12: rich truth (RegulonDB promoters and TF sites, terminators, REP, IS, rRNA, prophages).
- *Synechococcus elongatus* PCC 7942 or *Nostoc* PCC 7120: HIP1.
- *C. jejuni*: AT-rich.
- *B. cenocepacia*: GC-rich, multi-replicon.

**Nulls.**
- A genome simulated from each organism's own variable-order Markov model (order 5–8, both strands, same length).
- **The calibrated tree must make (almost) no splits below the composition levels on it.** This is the calibration
  test, and HASTUR's discipline applies unchanged.

**Readouts.**
1. For each annotation class, the best-matching node's F1 (precision and recall of member units against class
   units).
2. Leaf-partition normalised mutual information against each annotation layer.
3. The number of calibrated levels and nodes on the real genome vs the Markov null.
4. Runtime.

**Comparators.**
- k-means/SOM on 4-mer frequencies (BLSOM-like).
- Jensen–Shannon recursive segmentation (domains).
- Leiden on averaged gLM embeddings in 100-bp windows (a small, CPU-runnable model; Evo 2 SAE features if a GPU is
  available).
- The same tree using only contiguous k-mers. This ablation answers "do gapped patterns add anything?"

**Gate (decide before building v1).** The tree must:
- recover rRNA, REP/BIME, HIP1 and Chi-rich/strand classes, plus σ70 promoters and terminators as nodes, with F1
  clearly above the 4-mer baseline;
- stay quiet on the Markov null.

If promoters and terminators are not separated from generic intergenic sequence any better than by 4-mer
clustering, the "loose tendencies over tens of bases" premise does not hold at this resolution. Stop there.

---

## 8. (Superseded by Revision 1) Original organism plan

- **Novel organism, once the gate passes.**
  - A giant virus such as Mimivirus (1.2 Mb), with a known early promoter motif as a partial truth [bg].
  - Or a recently closed Asgard archaeon genome, such as *Prometheoarchaeum* (~4.5 Mb) [bg], where little non-coding
    biology is annotated.
  - The deliverable: a map, and a list of nodes that are calibrated, non-homologous, intergenic and unexplained by
    composition. These are candidate regulatory grammars for wet-lab follow-up.
- **Comparative.**
  - Build trees for related genomes, and match nodes by their pattern dictionaries.
  - A tendency that recurs across a genus is more likely functional.
  - The 21 local genomes are a pilot panel.
- **Pangenome or metagenome:** one tree over many genomes, so nodes become cross-species element classes.

---

## 9. Risks and open questions

1. **Power at the unit level.** A 100–200-bp unit holds few occurrences of any one gapped cell, so individual
   assignment is noisy. Node-level pooling carries the statistics, and soft assignment plus spatial HMM smoothing
   carries the per-unit labels. Unit width W is a registered parameter, with W = 64/128/256 compared on the gate.
2. **Mixed windows.** Several elements sit in one window. Hard trees will create "mixture" nodes, which is why §4.4
   is the planned v1.
3. **Coding dominance.** About 88% of a bacterial genome is coding, so the first splits will be about codon usage and
   frame. That is fine and expected; the intergenic subtree is the interesting part. A "non-coding only" mode is a
   simple input filter if wanted.
4. **Calibration subtleties.**
   - The null at each node must reproduce the node's context statistics *and* its unit-length distribution.
   - Overlapping windows must not be used for statistics.
   - Simulated-null artefacts are real artefacts (DEV_INSIGHTS §8.2): validate the generator.
5. **Meaning.** Calibration proves a node is not chance; it does not prove function. Annotation overlap,
   cross-species recurrence and positional regularity (fixed distance to start codons, for example) are the evidence
   of meaning, and each node report must show them.
6. **Prior art.** The quick search found nothing identical, but the space (composition maps, gLM SAEs, segmentation)
   is crowded. A focused search on "hierarchical clustering of genomic windows by k-mer/gapped-k-mer profiles" and on
   Evo 2 SAE prokaryotic feature atlases should be done before any claim of novelty.

---

## 10. Suggested first build (workstation only; the box is not needed)

1. **Rust:** per-unit gapped counting with streamed Σc and Σc² per cell, per node, extending `rust/prometheus`'s
   counting; a variable-order context model (BaMM back-off) and the expected-count computation.
2. **Python:** vocabulary selection, mixture split, null simulation and the calibration loop; then the tree driver
   and JSON output.
3. **Gate run:** *E. coli*, a cyanobacterium, and their Markov-null twins; readouts as in §7.
4. **Visualiser v0:** a static radial icicle plus edge bundling, with node panels, as one HTML page.

Estimated compute for the gate is minutes per genome. The work is mostly code, and per DEV_INSIGHTS §5.3 code time
is not the constraint.
