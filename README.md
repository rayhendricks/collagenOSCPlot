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
- **Add a gene** — type in the search box (symbol or WBGene ID) and pick from the list.
  Add multiple genes to overlay them; each gets its own color.
- **Remove** — click the `×` on a gene chip, or **Clear all**.
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
        ▼
        │  build_data.py  (join on stable WBGene ID, pack to JSON)
        ▼
data.json → data.js          loaded by index.html (Plotly dashboard)
```

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

# 5. build the dashboard data  ->  data.json, then wrap into data.js
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
| `build_data.py` | ✓ | Joins counts + annotation into `data.json` |
| `GSE130811_expr.tab.gz` | ✓ | Raw count matrix from GEO |
| `normalized_counts.tsv`, `data.json`, `*.gff.gz`, `gene_annotation.tsv` | — | Regenerable intermediates (git-ignored) |
