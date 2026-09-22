"""Small-multiples overview of every Teuscher 2019 collagen across the Meeuse 2020
time course (GSE130811), rendered as SVG -> PNG with headless Chrome.

    python3 plot_collagen_grid.py            # writes collagen_overview.html + .png
    python3 plot_collagen_grid.py --genes col-19,col-140,... --out late_collagens --cols 6
        # just those genes, in the given order, larger panels with hour ticks

Time is shown 20 °C-equivalent (25 °C hours x 1.5; see README "Time axis").
Each panel has its own linear y-scale (peak value printed top-right), so shapes are
comparable across genes but heights are not. Faint verticals = 20 °C molts
(Byerly et al. 1976: M1-M4 at 15, 23.5, 32.5, 45 h after hatching).
"""
import argparse, csv, json, os, subprocess, html, shutil

T20 = 1.5
MOLTS_20C = [15, 23.5, 32.5, 45]
ap = argparse.ArgumentParser()
ap.add_argument("--genes", help="comma-separated symbols (or a file, one per line); default = all collagens")
ap.add_argument("--out", default="collagen_overview", help="output basename (.html/.png)")
ap.add_argument("--cols", type=int, default=None)
ap.add_argument("--title", default=None)
ap.add_argument("--minimal", action="store_true", help="no title/legend/section headers/peak labels")
args = ap.parse_args()
subset = args.genes is not None
TICKS = subset                                # hour tick labels on larger panels
COLS = args.cols or (6 if subset else 12)
PW, PH, GAP = (300, 150, 10) if subset else (150, 84, 8)   # panel grid (CSS px)
PAD_L, PAD_R, PAD_T, PAD_B = (8, 8, 20, 20) if TICKS else (6, 6, 18, 6)   # plot area inside each panel
MARGIN, HEAD_H, SEC_H = (16, 0, 0) if args.minimal else (24, 92, 34)

C = dict(surface="#fcfcfb", panel="#ffffff", border="#e4e3de", grid="#ecebe6",
         text="#0b0b0b", text2="#52514e", muted="#8a8984",
         osc="#2a78d6", nonosc="#9a9993")    # osc = categorical slot 1; gray = de-emphasis

here = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(here, "data.js")).read()
D = json.loads(s[s.index("=") + 1:].rstrip().rstrip(";"))
ds = D["datasets"]["meeuse_5_48"]
hours = [h * T20 for h in ds["hours"]]
x0, x1 = hours[0], hours[-1]

rows = list(csv.DictReader(open(os.path.join(here, "collagens_teuscher2019.tsv")), delimiter="\t"))
genes = []
for r in rows:
    wb = r["wbgene"]
    if not wb or wb not in ds["series"]:
        continue
    v = ds["series"][wb]["v"]
    g = D["genes"][wb]
    genes.append(dict(sym=r["symbol"], cat=r["category"], v=v, vmax=max(v),
                      tmax=hours[v.index(max(v))], osc=g["osc"], phase=g["phase"]))

cut = [g for g in genes if g["cat"] == "Cuticular Collagens"]
if subset:
    raw = open(args.genes).read().split() if os.path.isfile(args.genes) else args.genes.split(",")
    want = []
    for n in (x.strip() for x in raw):
        n = n.replace("col-0", "col-")          # tolerate zero-padded names like col-08
        if n and n not in want:
            want.append(n)
    by = {g["sym"]: g for g in genes}
    missing = [n for n in want if n not in by]
    if missing:
        print("not found / not measured:", ", ".join(missing))
    picked = [by[n] for n in want if n in by]
    sections = [(f"Selected collagens, in the order given", picked)]
else:
  sections = [
    ("Cuticular collagens — oscillating (Meeuse 2020), ordered by peak phase",
     sorted([g for g in cut if g["osc"]], key=lambda g: g["phase"])),
    ("Cuticular collagens — not oscillating, ordered by time of peak expression",
     sorted([g for g in cut if not g["osc"]], key=lambda g: g["tmax"])),
    ("Non-cuticular collagens (basement membrane / other)",
     sorted([g for g in genes if g["cat"] == "Collagens"], key=lambda g: g["sym"])),
  ]

def fmt(n):
    return f"{n/1000:.0f}k" if n >= 10000 else (f"{n/1000:.1f}k" if n >= 1000 else f"{n:.0f}")

def panel(g, px, py):
    iw, ih = PW - PAD_L - PAD_R, PH - PAD_T - PAD_B
    sx = lambda t: px + PAD_L + (t - x0) / (x1 - x0) * iw
    top = g["vmax"] or 1
    sy = lambda y: py + PAD_T + ih - (y / top) * ih
    col = C["osc"] if g["osc"] else C["nonosc"]
    out = [f'<rect x="{px}" y="{py}" width="{PW}" height="{PH}" rx="4" fill="{C["panel"]}" stroke="{C["border"]}"/>']
    for m in MOLTS_20C:
        out.append(f'<line x1="{sx(m):.1f}" x2="{sx(m):.1f}" y1="{py+PAD_T}" y2="{py+PH-PAD_B}" stroke="{C["grid"]}" stroke-width="1"/>')
    out.append(f'<line x1="{px+PAD_L}" x2="{px+PW-PAD_R}" y1="{sy(0):.1f}" y2="{sy(0):.1f}" stroke="{C["border"]}"/>')
    pts = " ".join(f"{sx(t):.1f},{sy(y):.1f}" for t, y in zip(hours, g["v"]))
    out.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.5" stroke-linejoin="round" stroke-linecap="round"/>')
    if TICKS:
        for t in range(12, int(x1) + 1, 12):
            out.append(f'<text x="{sx(t):.1f}" y="{py+PH-6}" font-size="9" text-anchor="middle" fill="{C["muted"]}">{t}</text>')
    out.append(f'<text x="{px+6}" y="{py+13}" font-size="11" font-weight="600" fill="{C["text"]}">{html.escape(g["sym"])}</text>')
    right = fmt(g["vmax"]) + (f' · φ{g["phase"]:.0f}°' if g["osc"] and g["phase"] is not None else "")
    if not args.minimal:
        out.append(f'<text x="{px+PW-6}" y="{py+13}" font-size="9" text-anchor="end" fill="{C["text2"]}">{right}</text>')
    return "\n".join(out)

W = MARGIN * 2 + COLS * PW + (COLS - 1) * GAP
body, y = [], MARGIN + HEAD_H
for title, gs in sections:
    if not args.minimal:
        body.append(f'<text x="{MARGIN}" y="{y+20}" font-size="14" font-weight="600" fill="{C["text"]}">{html.escape(title)} '
                    f'<tspan fill="{C["text2"]}" font-weight="400">({len(gs)})</tspan></text>')
    y += SEC_H
    for i, g in enumerate(gs):
        r, c = divmod(i, COLS)
        body.append(panel(g, MARGIN + c * (PW + GAP), y + r * (PH + GAP)))
    y += ((len(gs) + COLS - 1) // COLS) * (PH + GAP) + (0 if args.minimal else 10)
H = y + MARGIN - (GAP if args.minimal else 0)

n = sum(len(gs) for _, gs in sections)
head = f'''
<text x="{MARGIN}" y="{MARGIN+18}" font-size="20" font-weight="700" fill="{C["text"]}">{html.escape(args.title or f"C. elegans collagens across larval development — {n} genes (Teuscher et al. 2019)")}</text>
<text x="{MARGIN}" y="{MARGIN+40}" font-size="12" fill="{C["text2"]}">mRNA-seq, Meeuse et al. 2020 (GSE130811), DESeq2-normalized counts · grown at 25 °C, time shown 20 °C-equivalent (25 °C h × 1.5): {x0:g}–{x1:g} h after plating · each panel has its own linear y-scale; top-right = peak count · phase{" · x-axis: 20 °C-equiv. hours" if TICKS else ""}</text>
<line x1="{MARGIN}" x2="{MARGIN+22}" y1="{MARGIN+62}" y2="{MARGIN+62}" stroke="{C["osc"]}" stroke-width="2.5"/>
<text x="{MARGIN+28}" y="{MARGIN+66}" font-size="12" fill="{C["text"]}">oscillating (Meeuse 2020)</text>
<line x1="{MARGIN+200}" x2="{MARGIN+222}" y1="{MARGIN+62}" y2="{MARGIN+62}" stroke="{C["nonosc"]}" stroke-width="2.5"/>
<text x="{MARGIN+228}" y="{MARGIN+66}" font-size="12" fill="{C["text"]}">not oscillating</text>
<line x1="{MARGIN+360}" x2="{MARGIN+360}" y1="{MARGIN+54}" y2="{MARGIN+70}" stroke="#c9c8c2" stroke-width="1.5"/>
<text x="{MARGIN+368}" y="{MARGIN+66}" font-size="12" fill="{C["text"]}">20 °C molts M1–M4 (15, 23.5, 32.5, 45 h; Byerly 1976)</text>
{"" if subset else f'<text x="{MARGIN+720}" y="{MARGIN+66}" font-size="12" fill="{C["text2"]}">Not measured: bli-6 (= col-112), col-182 (pseudogene)</text>'}
'''
if args.minimal:
    head = ""
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
       f'font-family="Inter, Helvetica, Arial, sans-serif"><rect width="{W}" height="{H}" fill="{C["surface"]}"/>'
       + head + "\n".join(body) + "</svg>")
page = f'<!doctype html><html><head><meta charset="utf-8"><title>Collagen overview</title>' \
       f'<style>html,body{{margin:0;background:{C["surface"]}}}</style></head><body>{svg}</body></html>'
hp = os.path.join(here, args.out + ".html"); open(hp, "w").write(page)

png = os.path.join(here, args.out + ".png")
tmp = os.environ.get("GRID_TMP", "/tmp")   # scratch dir for Chrome profile + screenshot
env = {k: v for k, v in os.environ.items() if k != "TMPDIR"}   # a long TMPDIR crashes Chrome
prof, shot = os.path.join(tmp, "collagen_grid_chrome"), os.path.join(tmp, args.out + ".png")
subprocess.run(["google-chrome", "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
                f"--user-data-dir={prof}", "--force-device-scale-factor=2",
                f"--window-size={W},{H}", f"--screenshot={shot}", "file://" + hp],
               check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
shutil.copyfile(shot, png)   # Chrome may be sandboxed from writing into the repo directly
print(f"{n} panels -> {png} ({W}x{H} CSS px @2x)")
