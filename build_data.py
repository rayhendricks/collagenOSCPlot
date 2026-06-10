#!/usr/bin/env python3
"""Join DESeq2-normalized counts with WBcel235 annotation into one compact JSON
for the dashboard."""
import json, re, csv

# --- annotation: wbgene -> (symbol, chrom, start, end, strand, biotype) ---
ann = {}
with open("gene_annotation.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        ann[row["wbgene"]] = (row["symbol"], row["chrom"], int(row["start"]),
                              int(row["end"]), row["strand"], row["biotype"])

# --- PRIMARY oscillation calls: Meeuse et al. 2020 Dataset EV1 ---
#     wbgene -> (is_oscillator, log2_amplitude, peak_phase_deg)
def _f(s):
    return float(s) if s not in ("", "nan", "NaN") else None

meeuse = {}
with open("meeuse_osc.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        meeuse[row["wbgene"]] = (row["osc"] == "1", _f(row["amp"]), _f(row["phase"]))

# --- cross-check: our computed oscillation calls (peaks/troughs heuristic) ---
osc_c = {}
with open("oscillation.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        osc_c[row["wbgene"]] = row["osc"] == "1"

# --- transcription factors (GO:0003700-defined): set of wbgene ---
tf = set()
with open("tf.tsv") as f:
    next(f)
    for line in f:
        tf.add(line.split("\t")[0])

# --- normalized counts ---
with open("normalized_counts.tsv") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]  # drop 'wbgene'

    main_idx, main_hours = [], []
    rep_idx, rep_hours = [], []
    for i, name in enumerate(sample_cols):
        m = re.match(r"^(\d+)hr(\.2)?$", name)
        hr = int(m.group(1))
        if m.group(2):  # replicate (.2)
            rep_idx.append(i); rep_hours.append(hr)
        else:
            main_idx.append(i); main_hours.append(hr)

    # genes  = dataset-independent annotation + osc/TF flags (computed once, shared)
    # series = per-dataset time course (v = main, r = replicates), keyed by WBGene ID
    genes = {}
    series = {}
    index = []
    for line in f:
        parts = line.rstrip("\n").split("\t")
        wb = parts[0]
        vals = [float(x) for x in parts[1:]]
        sym, chrom, start, end, strand, bt = ann.get(
            wb, (wb, "?", 0, 0, ".", "unknown"))
        m = meeuse.get(wb)                       # Meeuse 2020 (primary)
        is_osc = m[0] if m else False
        amp = m[1] if m else None
        phase = m[2] if m else None
        is_oscC = osc_c.get(wb, False)           # our computed (cross-check)
        is_tf = wb in tf
        genes[wb] = {
            "sym": sym, "chrom": chrom, "start": start, "end": end,
            "strand": strand, "bt": bt,
            "osc": is_osc, "amp": amp, "phase": phase, "oscC": is_oscC, "tf": is_tf,
        }
        series[wb] = {
            "v": [round(vals[i], 2) for i in main_idx],
            "r": [round(vals[i], 2) for i in rep_idx],
        }
        index.append({"wb": wb, "sym": sym,
                      "o": 1 if is_osc else 0, "t": 1 if is_tf else 0})

# --- second design: ribosome footprinting (Hendriks et al. 2014, GSE52905) ---
#     Continuous development, N2, 18-36 h, every 2 h, 10 timepoints (no replicates).
#     The supplementary matrix is log2 depth-normalized footprint counts; we store
#     2**x so the values are linear footprint abundance and the dashboard's
#     linear/log toggle behaves exactly as it does for the mRNA design. This is a
#     SEPARATE time base (continuous-development hours) and a different quantity
#     (translation, not transcript) — never overlaid on the mRNA axis.
import gzip
fp_series = {}
fp_hours = []
with gzip.open("GSE52905_footprint_normalized.txt.gz", "rt") as f:
    fp_cols = f.readline().rstrip("\n").split("\t")[1:]
    fp_hours = [int(re.search(r"_(\d+)h", c).group(1)) for c in fp_cols]
    for line in f:
        p = line.rstrip("\n").split("\t")
        wb = p[0]
        if wb not in genes:        # keep searchable genes only (shared annotation)
            continue
        fp_series[wb] = {"v": [round(2.0 ** float(x), 2) for x in p[1:]]}
print(f"footprint genes mapped to annotation: {len(fp_series)}  hours: {fp_hours}")

out = {
    "meta": {
        "series": "GSE130811",
        "title": "C. elegans extended wild-type (N2) developmental time course, 5-48h @25C",
        "normalization": "DESeq2 median-of-ratios (size-factor) normalized counts",
        "assembly": "WBcel235 / ce11 (RefSeq GCF_000002985.6)",
        "n_genes": len(genes),
        "n_osc": sum(1 for g in genes.values() if g["osc"]),
        "n_oscC": sum(1 for g in genes.values() if g["oscC"]),
        "n_tf": sum(1 for g in genes.values() if g["tf"]),
        "osc_source": "Meeuse et al. 2020 (Mol Syst Biol), Dataset EV1 — cosine fit",
    },
    "genes": genes,
    "index": index,
    # selectable designs; each bundle owns its own axes/units so incomparable
    # time courses (e.g. a future Ribo-seq dauer-exit series) never share an axis.
    "default_dataset": "meeuse_5_48",
    "datasets": {
        "meeuse_5_48": {
            "meta": {
                "label": "Continuous development — mRNA (Meeuse 2020)",
                "series": "GSE130811",
                "assay": "mRNA-seq",
                "units": "norm. counts",
                "timeref": "hours after plating (25°C)",
            },
            "hours": main_hours,
            "rep_hours": rep_hours,
            "series": series,
        },
        "footprint_contDev": {
            "meta": {
                "label": "Continuous development — ribosome footprinting (Hendriks 2014)",
                "series": "GSE52905",
                "assay": "Ribo-seq (ribosome footprinting)",
                "units": "norm. footprint",
                "timeref": "hours of continuous development (25°C)",
            },
            "hours": fp_hours,
            "rep_hours": [],
            "series": fp_series,
        },
    },
}
with open("data.json", "w") as f:
    json.dump(out, f, separators=(",", ":"))

print(f"genes: {len(genes)}  main timepoints: {main_hours}")
print(f"replicate timepoints: {rep_hours}")
n_unmapped = sum(1 for wb in genes if genes[wb]['chrom'] == '?')
print(f"unmapped to WBcel235 annotation: {n_unmapped}")
import os
print(f"data.json size: {os.path.getsize('data.json')/1e6:.1f} MB")
