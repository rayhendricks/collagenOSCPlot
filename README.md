# C. elegans developmental oscillator — expression dashboard

## **[LAUNCH THE LIVE DASHBOARD](https://rayhendricks.github.io/collagenOSCPlot/)**

### https://rayhendricks.github.io/collagenOSCPlot/

[![Open Dashboard](https://img.shields.io/badge/OPEN_LIVE_DASHBOARD-rayhendricks.github.io%2FcollagenOSCPlot-2ea44f?style=for-the-badge)](https://rayhendricks.github.io/collagenOSCPlot/)

*Runs in any browser — no install, no account. Search a gene, overlay several, toggle linear/log.*

---

An interactive browser dashboard for exploring **DESeq2-normalized RNA-seq expression**
across the *C. elegans* extended wild-type (N2) developmental time course
(**GEO [GSE130811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130811)**,
a sub-series of the SuperSeries
[GSE133576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133576),
*"State transitions of a developmental oscillator"*).

Search any gene by symbol (`lin-42`, `col-19`) or WBGene ID, overlay several genes at once,
and toggle the y-axis between **linear** (absolute amplitude) and **log₁₀** (fold-change rhythm).

![timepoints](https://img.shields.io/badge/timepoints-5–48h%20hourly-blue) ![genes](https://img.shields.io/badge/genes-20%2C392-blue) ![normalization-DESeq2-blue](https://img.shields.io/badge/normalization-DESeq2-blue)

---

## Quick start (just view it)

No build step required — the normalized data is bundled in `data.js`:

```bash
open index.html        # macOS
# or: python -m http.server  then visit http://localhost:8000
```

Opening the file directly works because the data is loaded via `<script src="data.js">`
(not `fetch`), so there are no `file://` CORS issues.

### Using the dashboard
- **Dataset** — pick the experimental design from the dropdown. Two are wired up:
  the **mRNA** continuous-development time course (Meeuse 2020, 5–48 h) and the
  **ribosome-footprinting** continuous-development course (Hendriks 2014, 18–36 h).
  Switching reswaps the axes, units, and timepoints and redraws your selected genes.
  Each design is **normalized independently** and on a **different time base**, so
  absolute heights are not comparable across datasets — the rhythm and phase are.
  A gene not measured in the active dataset shows *"no data in this dataset"* on its chip.
- **Add a gene** — type in the search box (symbol or WBGene ID) and pick from the list.
  Add multiple genes to overlay them; each gets its own color.
- **Remove** — click the `×` on a gene chip, or **Clear all**.
- **Oscillators only** — restrict the search to the 3,739 genes classified as rhythmic by
  Meeuse et al. 2020 (see *Oscillator calls* below). Rhythmic genes carry an `osc φN°`
  badge (N = peak phase in degrees), and their amplitude + phase show on the gene chip.
- **Transcription factors only** — restrict the search to genes with DNA-binding
  transcription-factor activity (see *Transcription factors* below); these carry a
  `TF` badge. Combine with **oscillators only** to find the rhythmic TFs (the
  oscillator's own regulators — 125 genes).
- **Replicate points** — 7 timepoints (37, 38, 39, 41, 45, 47, 48 h) have a second
  replicate, shown as open diamonds. Toggle with the checkbox.
- **Linear vs Log₁₀** — linear shows how big a peak is; log shows the rhythm and is the
  honest view for sharp oscillators like `lin-42`.

---

## The dataset

| | |
|---|---|
| Series | GSE130811 (sub-series of GSE133576) |
| Organism | *Caenorhabditis elegans*, strain N2 (wild type) |
| Assay | mRNA-seq, Illumina HiSeq 2500 |
| Design | Single developmental time course, **5–48 h** post-plating at 25 °C, sampled hourly (44 timepoints) + 7 replicate samples = **51 samples** |
| Genes | 20,392 (protein-coding-centric matrix) |

> **Scope:** this is a **visualization** tool. The series is one trajectory with ~n=1 per
> timepoint, so it supports describing each gene's expression over time — not statistical
> differential-expression or rhythmicity testing (which would need replication or a periodic model).

---

## Pipeline

```
GSE130811_expr.tab.gz        raw integer counts (20,392 genes × 51 samples) from GEO
        │  normalize_deseq2.R (DESeq2 median-of-ratios size factors)
        ▼
normalized_counts.tsv  +  size_factors.tsv
        │
WBcel235_genomic.gff.gz      latest RefSeq assembly (GCF_000002985.6)
        │  awk → gene_annotation.tsv   (WBGene ID ↔ symbol ↔ coords ↔ biotype)
        │
normalized_counts.tsv
        │  detect_oscillators.py → oscillation.tsv      (computed osc cross-check)
meeuse_EV1.xlsx (Dataset EV1)
        │  parse_meeuse.py      → meeuse_osc.tsv        (PRIMARY osc + amplitude + phase)
wb.gaf.gz + go-basic.obo
        │  detect_tfs.py        → tf.tsv                (transcription factors, GO:0003700)
        │
GSE52905_..._Footprint_normalized.txt.gz   ribosome footprinting (Hendriks 2014)
        │  (read directly by build_data.py as a second selectable dataset)
        ▼
        │  build_data.py  (shared gene annotation + osc + TF; one time-course
        │                  bundle per dataset, joined on WBGene ID)
        ▼
data.json → data.js          loaded by index.html (Plotly dashboard)
```

### Datasets (the dashboard's dropdown)
The dashboard hosts several **independently-normalized designs** rather than merging
incomparable time courses. A gene's annotation and osc/TF flags are stored once
(dataset-independent); each dataset carries only its own time series, keyed on WBGene ID.

| Dataset | Source | Assay | Window | Quantity |
|---|---|---|---|---|
| Continuous development — mRNA | GSE130811 (Meeuse 2020) | mRNA-seq | 5–48 h post-plating, hourly | DESeq2 norm. counts |
| Continuous development — footprinting | GSE52905 (Hendriks 2014) | Ribo-seq | 18–36 h cont. dev., every 2 h | log2 depth-norm. footprint counts (stored linearized) |

Both are N2 @25 °C, but the **time bases differ** (hours after plating vs. hours of
continuous development) and the **measured quantity differs** (transcript abundance vs.
ribosome occupancy), so they are never overlaid on one axis — you switch between them.

### Oscillator calls — primary: Meeuse et al. 2020 (`parse_meeuse.py`)
The dashboard's `osc` flag is the **authoritative, peer-reviewed classification** from
[Meeuse et al. 2020](https://doi.org/10.15252/msb.20209498) — the paper this dataset comes
from — supplied in their **Dataset EV1** (`meeuse_EV1.xlsx`, included; CC-BY 4.0). They
classified **3,739 genes (24% of expressed)** as oscillating by **cosine fitting** (fixed
7-h period, on t = 10–25 h where the period is most stable), with cut-offs **log₂
amplitude ≥ 0.5 (≥ 2-fold)** and **P ≤ 0.01**. Crucially, every gene also carries a
**log₂ amplitude** and a **peak phase** (degrees), so the dashboard shows *when* in the
cycle each gene peaks (the `osc φN°` badge) — not just whether it oscillates. `parse_meeuse.py`
flattens the table to `meeuse_osc.tsv`, joined on stable WBGene ID (3,739/3,739 map to our
matrix).

### Oscillator calls — cross-check: computed (`detect_oscillators.py`)
We also **compute oscillation from the time course itself** (stored as `oscC`), as a
reproducible, dependency-free cross-check. The target signal is the *C. elegans*
**molt-cycle oscillation**: genes that rise and fall rhythmically with a period of roughly
one larval stage (~7–8 h). Three properties of this real signal break the textbook
periodicity tests, and shaped the method:

| Property of the real signal | Why naive methods fail |
|---|---|
| Peaks are **sharp/spiky** (non-sinusoidal) | A single FFT bin sees the energy spread into harmonics → `lin-42` scores ~0.05 |
| The period **drifts** (lengthens: 7 → 7 → 8 h) | Fixed-lag autocorrelation can't align all cycles → weak/negative score |
| Only **~4–6 cycles** fit in 5–48 h | Too few cycles for stable spectral/autocorrelation statistics |

So instead of scoring periodicity in the frequency domain, the detector finds the
biology directly — *regularly spaced peaks that return to baseline.* For each gene
(operating on the 44 single-replicate hourly samples, 5–48 h):

**1. Log transform.** `y = log2(count + 1)`. Expression amplitude is multiplicative, so
log space makes a fold-change rhythm look like a constant-amplitude wave; the `+1`
keeps zeros finite.

**2. Detrend with a *local* moving average.** Subtract a centered 7-h rolling mean
(`window ≈ one period`) to remove the slow developmental drift while preserving the
~7-h cycle: `resid = y − rolling_mean(y, 7)`. The choice of a **local** smoother is the
key trick — a global polynomial fit *rings* (Gibbs-style wiggles) when it tries to follow
a one-time step, which fakes an oscillation and made switch genes like `col-19` score as
cyclers. A local average has nothing to ring with: on a step it just tracks the step.

**3. Find prominent peaks _and_ troughs.** On `resid`, a point is a peak if it is a local
maximum rising at least **30 % of the gene's amplitude** (`amp = max(resid) − min(resid)`)
above the floor; peaks closer than **4 h** are merged (keep the taller). Troughs are the
same on `−resid`. Requiring **both** peaks and troughs is what separates an oscillator
(repeatedly returns to baseline: up–down–up–down) from a **switch** like `col-19` (rises
once and stays high — many peaks possible, but it never comes back down repeatedly).

**4. Estimate the period.** `period = median of the gaps (in hours) between consecutive
peaks` (the molt-cycle period; the dashboard now shows Meeuse's peak *phase* instead).

**5. Call it.** A gene is flagged **oscillating** when *all* hold:

| Criterion | Threshold | Purpose |
|---|---|---|
| max normalized count | ≥ 16 | ignore genes that never rise above noise |
| amplitude `amp` | ≥ 2.5 log₂ (~5.7-fold) | require a real, large swing |
| number of peaks | ≥ 3 | several cycles, not one bump |
| number of troughs | ≥ 3 | repeatedly returns to baseline (excludes switches) |
| median peak interval (period) | 6–10 h | the molt-cycle band |

**Calibration & validation.** These thresholds yield **3,062** computed oscillators. On
spot-checks the method recovers **8/8** canonical cyclers (`lin-42, mlt-10, dpy-13, bli-1,
qua-1, mlt-9, noah-1, grd-3`) and rejects **4/4** controls (`act-1` — too low amplitude;
`col-19` — a switch, period out of band; `tba-1`, `ama-1` — flat).

**Agreement with Meeuse.** Of the 3,062 computed calls, **2,071 (68%) are also Meeuse
oscillators**; conversely the computed set recovers **55%** of Meeuse's 3,739. The gap is
expected and instructive: (a) our amplitude cut-off (≥ 5.7-fold) is **far stricter** than
Meeuse's (≥ 2-fold), so we drop ~1,668 genuine low-amplitude cyclers (e.g. `blmp-1`, a key
oscillator-regulating TF, amplitude 0.68); (b) our peak/trough method catches sharp,
**non-sinusoidal** oscillators that Meeuse's cosine fit can miss, contributing some of the
991 computed-only calls. The computed flag is a deliberately **high-confidence,
conservative** heuristic, not a statistical rhythmicity test — which is exactly why the
**Meeuse list is the primary `osc` flag** and the computed one is kept only as a
cross-check.

### Transcription factors (`detect_tfs.py`)
TF status is **not** inferred from expression — it comes from curated annotation. A gene
is flagged a TF if its product has **DNA-binding transcription factor activity**
(Gene Ontology **GO:0003700**), i.e. it binds DNA sequence-specifically to regulate
transcription. The script parses `go-basic.obo` to collect GO:0003700 **and all its
descendant terms** (the RNA-pol-II-specific children etc.), then selects every
*C. elegans* gene with an `enables` annotation to any of them in WormBase's GO file
(`wb.gaf.gz`).

**Result:** **632** TFs in the matrix. Validates 12/13 spot-checked TFs (`daf-16, pha-4,
skn-1, blmp-1, nhr-23/25, ...`; the lone miss, `lin-14`, genuinely lacks a GO DNA-binding
annotation) and 0/7 non-TFs. This is the **curated, high-confidence** definition — more
conservative than the ~900 *computationally predicted* TFs in compendia like wTF, by
design. Source: GO consortium GAF + go-basic ontology.

### Normalization detail
`normalize_deseq2.R` drops the gene-length `width` column, builds a `DESeqDataSet`
with `design = ~1` (intercept only — DESeq2 is used **purely as a normalization engine**,
no GLM/testing), runs `estimateSizeFactors()`, and exports `counts(dds, normalized = TRUE)`.
Median-of-ratios rescales each sample by a single size factor (observed range **0.54–1.33**),
making samples comparable without distorting any gene's temporal shape.

### Annotation detail
Gene symbols and coordinates come from the **latest RefSeq assembly, WBcel235 / ce11**
(`GCF_000002985.6`). The join key is the **stable WBGene ID**, so it is reliable even though
the count matrix was originally built on ce10. 1,208 genes (retired/merged WBGene IDs) have
no current WBcel235 coordinate but remain searchable/plottable by ID.

---

## Reproducing from scratch

Requires [mamba](https://github.com/conda-forge/miniforge).

```bash
# 1. environment
mamba create -n collagen -c conda-forge -c bioconda -y \
  r-base bioconductor-deseq2 r-data.table python=3.11 pandas

# 2. raw counts (already included as GSE130811_expr.tab.gz)
curl -sL "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130811&format=file&file=GSE130811%5Fexpr%5FmRNA%5FCE10%5Fcoding%2Etab%2Egz" \
  -o GSE130811_expr.tab.gz

# 3. DESeq2 normalization  ->  normalized_counts.tsv, size_factors.tsv
mamba run -n collagen Rscript normalize_deseq2.R

# 4. RefSeq annotation  ->  gene_annotation.tsv
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.gff.gz" \
  -o WBcel235_genomic.gff.gz
#   (parsing awk one-liner is documented in build steps / commit history)

# 5. oscillator calls  ->  oscillation.tsv (computed cross-check) + meeuse_osc.tsv (primary)
mamba run -n collagen python detect_oscillators.py
#   meeuse_EV1.xlsx is included; if missing, download Dataset EV1 from the paper
#   (https://doi.org/10.15252/msb.20209498, Supplementary Information) and save as meeuse_EV1.xlsx
mamba run -n collagen python parse_meeuse.py

# 6. transcription-factor calls  ->  tf.tsv
curl -sL "http://current.geneontology.org/annotations/wb.gaf.gz" -o wb.gaf.gz
curl -sL "http://purl.obolibrary.org/obo/go/go-basic.obo" -o go-basic.obo
mamba run -n collagen python detect_tfs.py

# 7. second dataset: ribosome footprinting (Hendriks 2014, GSE52905)
curl -sL "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE52905&format=file&file=GSE52905%5Fce%5FgeneExpression%5FFootprint%5Fnormalized%2Etxt%2Egz" \
  -o GSE52905_footprint_normalized.txt.gz

# 8. build the dashboard data  ->  data.json, then wrap into data.js
mamba run -n collagen python build_data.py
{ printf 'window.DATA='; cat data.json; } > data.js
```

---

## Files

| File | Tracked | Purpose |
|---|:--:|---|
| `index.html` | ✓ | The dashboard (Plotly + vanilla JS) |
| `data.js` | ✓ | Bundled normalized expression + annotation (`window.DATA`) |
| `normalize_deseq2.R` | ✓ | DESeq2 size-factor normalization |
| `parse_meeuse.py` | ✓ | Flattens Meeuse Dataset EV1 → primary osc + amplitude + phase |
| `meeuse_EV1.xlsx` | ✓ | Meeuse et al. 2020 Dataset EV1 (CC-BY 4.0) |
| `detect_oscillators.py` | ✓ | Computes osc cross-check from the time course |
| `detect_tfs.py` | ✓ | Flags transcription factors from GO:0003700 |
| `build_data.py` | ✓ | Builds shared annotation + one time-course bundle per dataset into `data.json` |
| `GSE130811_expr.tab.gz` | ✓ | Raw count matrix from GEO |
| `normalized_counts.tsv`, `oscillation.tsv`, `meeuse_osc.tsv`, `tf.tsv`, `data.json`, `*.gff.gz`, `*.gaf.gz`, `go-basic.obo`, `GSE52905_footprint_*.txt.gz`, `gene_annotation.tsv` | — | Regenerable intermediates (git-ignored) |
