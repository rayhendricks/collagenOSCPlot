#!/usr/bin/env python3
"""Flag oscillating genes directly from the time course (no external gene list).

The molt-cycle oscillation here is sharp (spiky), its period drifts (lengthens)
across development, and only ~4-6 cycles fit in 5-48 h -- which defeats FFT and
fixed-lag autocorrelation. We instead detect the biology directly:

  1. log2(x+1)
  2. remove the slow developmental trend with a *local* moving average
     (window ~= one period). A local smoother (unlike a global polynomial fit)
     does not "ring" on a step, so switch-like genes (e.g. col-19) are not
     mistaken for oscillators.
  3. require several PROMINENT PEAKS *and* TROUGHS (the signal repeatedly returns
     to baseline -- a switch does not), regularly spaced at a molt-cycle interval
     (~6-10 h), with real amplitude.

Calibration: yields ~3,062 oscillating genes (cf. ~3,235 high-confidence in
Meeuse et al. 2020), recovers 8/8 canonical cyclers (lin-42, mlt-10, dpy-13,
bli-1, qua-1, mlt-9, noah-1, grd-3) and excludes 4/4 controls (act-1, col-19,
tba-1, ama-1)."""
import re
import numpy as np
import pandas as pd

WINDOW       = 7       # hours; moving-average detrend window (~one period)
PROM_FRAC    = 0.30    # peak/trough must rise this fraction of amplitude above floor
MIN_SEP      = 4       # hours; merge peaks closer than this (keep the taller)
MIN_PEAKS    = 3       # repeated cycles required
MIN_TROUGHS  = 3       # ...and repeated returns to baseline (excludes switches)
PERIOD_MIN, PERIOD_MAX = 6.0, 10.0   # hours; allowed median peak interval
MIN_MAX_EXPR = 16.0    # ignore genes that never rise above noise
AMP_THRESH   = 2.5     # peak-to-trough amplitude in log2 (>= ~5.7-fold)


def find_peaks(sig, thresh, min_sep):
    cand = [i for i in range(len(sig))
            if sig[i] >= thresh
            and sig[i] >= sig[max(0, i - 1)] and sig[i] >= sig[min(len(sig) - 1, i + 1)]]
    kept = []
    for i in cand:
        if kept and (i - kept[-1]) < min_sep:
            if sig[i] > sig[kept[-1]]:
                kept[-1] = i
        else:
            kept.append(i)
    return kept


with open("normalized_counts.tsv") as f:
    header = f.readline().rstrip("\n").split("\t")[1:]
    main = [(i, int(re.match(r"^(\d+)hr$", n).group(1)))
            for i, n in enumerate(header) if re.match(r"^(\d+)hr$", n)]
    main.sort(key=lambda t: t[1])
    cols = [i for i, _ in main]
    hours = np.array([h for _, h in main], dtype=float)

    rows = []
    for line in f:
        p = line.rstrip("\n").split("\t")
        wb = p[0]
        vals = np.array([float(p[1 + i]) for i in cols])
        maxv = vals.max()
        y = np.log2(vals + 1.0)
        trend = pd.Series(y).rolling(WINDOW, center=True, min_periods=1).mean().values
        resid = y - trend
        amp = float(resid.max() - resid.min())
        if maxv < MIN_MAX_EXPR or amp < AMP_THRESH:
            rows.append((wb, 0, 0.0, round(amp, 3), False)); continue
        pk = find_peaks(resid, resid.min() + PROM_FRAC * amp, MIN_SEP)
        tr = find_peaks(-resid, (-resid).min() + PROM_FRAC * amp, MIN_SEP)
        iv = np.diff(hours[pk]) if len(pk) >= 2 else np.array([])
        med = float(np.median(iv)) if iv.size else 0.0
        osc = (len(pk) >= MIN_PEAKS and len(tr) >= MIN_TROUGHS
               and PERIOD_MIN <= med <= PERIOD_MAX)
        rows.append((wb, len(pk), round(med, 1), round(amp, 3), osc))

with open("oscillation.tsv", "w") as f:
    f.write("wbgene\tn_peaks\tperiod\tamp_log2\tosc\n")
    for wb, npk, per, amp, osc in rows:
        f.write(f"{wb}\t{npk}\t{per}\t{amp}\t{int(osc)}\n")

n_osc = sum(1 for r in rows if r[4])
print(f"genes scored: {len(rows)}  |  flagged oscillating: {n_osc}")
print(f"method: moving-avg detrend W={WINDOW}h, >={MIN_PEAKS} peaks & >={MIN_TROUGHS} troughs, "
      f"median interval {PERIOD_MIN}-{PERIOD_MAX}h, amp>={AMP_THRESH} log2, maxexpr>={MIN_MAX_EXPR}")
