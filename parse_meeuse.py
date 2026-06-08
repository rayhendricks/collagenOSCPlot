#!/usr/bin/env python3
"""Convert Meeuse et al. 2020 Dataset EV1 (the authoritative oscillating-gene
classification for this dataset) into a flat TSV keyed by WBGene ID.

Source: Mol Syst Biol 16:e9498, Dataset EV1 (manually downloaded as
meeuse_EV1.xlsx). 3,739 genes classified "Osc" by cosine fitting (fixed 7-h
period on t=10-25 h, log2 amplitude >= 0.5 i.e. >=2-fold, P <= 0.01). Each gene
carries a log2 amplitude and a peak phase (degrees)."""
import pandas as pd

df = pd.read_excel("meeuse_EV1.xlsx", sheet_name=0)
out = df[["WB_ID", "Class", "OscAmplitude", "PeakPhase"]].copy()
out["osc"] = (out["Class"] == "Osc").astype(int)
out["amp"] = out["OscAmplitude"].round(3)
out["phase"] = out["PeakPhase"].round(1)
out[["WB_ID", "osc", "amp", "phase"]].rename(columns={"WB_ID": "wbgene"}).to_csv(
    "meeuse_osc.tsv", sep="\t", index=False)

print(f"rows: {len(out)} | Osc: {int(out['osc'].sum())} | nonOsc: {int((1-out['osc']).sum())}")
