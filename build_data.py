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

# --- oscillation calls: wbgene -> (is_oscillator, period_hours) ---
osc = {}
with open("oscillation.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        osc[row["wbgene"]] = (row["osc"] == "1", float(row["period"]))

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

    genes = {}
    index = []
    for line in f:
        parts = line.rstrip("\n").split("\t")
        wb = parts[0]
        vals = [float(x) for x in parts[1:]]
        sym, chrom, start, end, strand, bt = ann.get(
            wb, (wb, "?", 0, 0, ".", "unknown"))
        is_osc, period = osc.get(wb, (False, 0.0))
        is_tf = wb in tf
        genes[wb] = {
            "sym": sym, "chrom": chrom, "start": start, "end": end,
            "strand": strand, "bt": bt,
            "osc": is_osc, "period": period, "tf": is_tf,
            "v": [round(vals[i], 2) for i in main_idx],
            "r": [round(vals[i], 2) for i in rep_idx],
        }
        index.append({"wb": wb, "sym": sym,
                      "o": 1 if is_osc else 0, "t": 1 if is_tf else 0})

out = {
    "meta": {
        "series": "GSE130811",
        "title": "C. elegans extended wild-type (N2) developmental time course, 5-48h @25C",
        "normalization": "DESeq2 median-of-ratios (size-factor) normalized counts",
        "assembly": "WBcel235 / ce11 (RefSeq GCF_000002985.6)",
        "n_genes": len(genes),
        "n_osc": sum(1 for g in genes.values() if g["osc"]),
        "n_tf": sum(1 for g in genes.values() if g["tf"]),
    },
    "hours": main_hours,
    "rep_hours": rep_hours,
    "genes": genes,
    "index": index,
}
with open("data.json", "w") as f:
    json.dump(out, f, separators=(",", ":"))

print(f"genes: {len(genes)}  main timepoints: {main_hours}")
print(f"replicate timepoints: {rep_hours}")
n_unmapped = sum(1 for wb in genes if genes[wb]['chrom'] == '?')
print(f"unmapped to WBcel235 annotation: {n_unmapped}")
import os
print(f"data.json size: {os.path.getsize('data.json')/1e6:.1f} MB")
