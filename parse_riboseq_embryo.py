#!/usr/bin/env python3
"""Extract the wild-type (N2) early-embryo time courses from GSE281412 (Cenik lab,
"Low input ribosome profiling … during early C. elegans embryogenesis") into two
gene x sample TSVs the dashboard can load:

    embryo_rnaseq.tsv    total RNA-seq        (whole-transcript counts)
    embryo_riboseq.tsv   ribosome profiling   (CDS footprint counts, RiboITP)

The processed data are RiboPy `.ribo` containers (HDF5). Two series-level files:
  * GSE281412.HDF5    — the RNA-seq batches; each experiment carries an `rnaseq`
                        quantification (n_transcripts x 5 transcript regions).
  * GSE281412_1.HDF5  — the ribosome-profiling experiments; `region_counts` is
                        (n_lengths * n_transcripts, 5 regions), read-length-major.

We keep only the WT_* experiments (the OMA-1 gain-of-function mutant is ignored),
sum to whole-transcript (RNA) or CDS (Ribo) counts, collapse isoforms to the stable
WBGene ID, and depth-normalize each sample to counts-per-million so samples on the
tiny early-embryo inputs are comparable. Columns are named c{cells}_{rep}
(cells in {1,2,4,8}); the dashboard build groups them into per-stage means + reps.
"""
import re, csv, urllib.request, os
from collections import defaultdict
import numpy as np
import h5py

RNA_FILE  = "GSE281412.HDF5"
RIBO_FILE = "GSE281412_1.HDF5"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281412/suppl/{}")

# RiboPy region order: [UTR5, UTR5_junction, CDS, UTR3_junction, UTR3]
CDS = 2
N_LENGTHS = 20          # read lengths 21..40 nt (length_min/max in the file attrs)

_gene_re = re.compile(r"gene:(WBGene\d+)")


def _download(fn):
    if not os.path.exists(fn):
        print(f"  downloading {fn} …")
        urllib.request.urlretrieve(URL.format(fn), fn)


def _gene_ids(h):
    """Transcript index -> stable WBGene ID (isoform suffix stripped)."""
    out = []
    for raw in h["reference/reference_names"][:]:
        m = _gene_re.search(raw.decode())
        out.append(m.group(1) if m else None)
    return out


def _stage(name):
    """WT_1cell_B11_1 / WT_1-cell_1 -> cell count (1,2,4,8)."""
    return int(re.search(r"_(\d+)-?cell", name).group(1))


def _collapse_to_genes(per_transcript, gene_ids):
    """Sum transcript-level counts to gene level (skip transcripts w/o WBGene)."""
    g = defaultdict(float)
    for v, wb in zip(per_transcript, gene_ids):
        if wb is not None:
            g[wb] += float(v)
    return g


def _write(path, gene_cols, col_names):
    genes = sorted(set().union(*[c.keys() for c in gene_cols.values()]))
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["wbgene"] + col_names)
        for wb in genes:
            w.writerow([wb] + [f"{gene_cols[c].get(wb, 0.0):.2f}" for c in col_names])
    print(f"  wrote {path}: {len(genes)} genes x {len(col_names)} samples")


def _cpm(gene_counts):
    """counts-per-million within one sample (depth normalization)."""
    tot = sum(gene_counts.values()) or 1.0
    return {wb: v * 1e6 / tot for wb, v in gene_counts.items()}


def _ordered_cols(stage_of):
    """Stage-then-replicate column names: c1_1, c1_2, …, c8_n."""
    by_stage = defaultdict(list)
    for exp, st in stage_of.items():
        by_stage[st].append(exp)
    cols, mapping = [], {}
    for st in sorted(by_stage):                       # 1,2,4,8
        for i, exp in enumerate(sorted(by_stage[st]), 1):
            name = f"c{st}_{i}"
            cols.append(name)
            mapping[name] = exp
        print(f"  {st}-cell: {len(by_stage[st])} reps")
    return cols, mapping


def parse_rnaseq():
    print(RNA_FILE)
    _download(RNA_FILE)
    with h5py.File(RNA_FILE, "r") as h:
        gene_ids = _gene_ids(h)
        wt = {e: _stage(e) for e in h["experiments"] if e.startswith("WT")}
        cols, mapping = _ordered_cols(wt)
        gene_cols = {}
        for name, exp in mapping.items():
            rna = h[f"experiments/{exp}/rnaseq/rnaseq"][:]   # (n_transcripts, 5)
            per_tx = rna.sum(axis=1)                          # whole transcript
            gene_cols[name] = _cpm(_collapse_to_genes(per_tx, gene_ids))
    _write("embryo_rnaseq.tsv", gene_cols, cols)


def parse_riboseq():
    print(RIBO_FILE)
    _download(RIBO_FILE)
    with h5py.File(RIBO_FILE, "r") as h:
        gene_ids = _gene_ids(h)
        n_tx = len(gene_ids)
        wt = {e: _stage(e) for e in h["experiments"] if e.startswith("WT")}
        cols, mapping = _ordered_cols(wt)
        gene_cols = {}
        for name, exp in mapping.items():
            rc = h[f"experiments/{exp}/region_counts/region_counts"][:]
            # read-length-major: (n_lengths, n_transcripts, 5); sum over lengths,
            # keep the CDS region — the standard ribosome-footprint quantification.
            per_tx = rc.reshape(N_LENGTHS, n_tx, 5)[:, :, CDS].sum(axis=0)
            gene_cols[name] = _cpm(_collapse_to_genes(per_tx, gene_ids))
    _write("embryo_riboseq.tsv", gene_cols, cols)


if __name__ == "__main__":
    parse_rnaseq()
    parse_riboseq()
