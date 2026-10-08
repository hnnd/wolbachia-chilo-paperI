#!/usr/bin/env python3
"""Figure S6 — MAG quality of the two recovered Wolbachia genomes.

Panel a: CheckM2 completeness–contamination plane with the MIMAG contamination
         ceiling and the two removal steps that produced the deposited wCchiA set.
Panel b: best-hit amino-acid identity of the 89 predicted proteins on all removed
         sequences to the deposited six-contig wCchiA assembly.

Inputs (data/tables/):
  figI2_genome_quality.csv           -- per-assembly CheckM2 / assembly statistics
  figI2_removed_protein_identity.csv -- panel b bin counts

Panel a coordinates come from CheckM2 v1.1.0 runs on the exact FASTA sets named in
figI2_genome_quality.csv.  Panel b was computed as follows: proteins were predicted with
`prodigal -p meta` on each of the four sequences removed before deposition
(wCchiA_contig_04, wCchiA_contig_08, NODE_12_length_17695_cov_9.432853_RagTag,
ptg000003l_RagTag; the `-p meta` flag is required because single mode needs >= 20 kb),
then searched with `diamond blastp` against the 1,326 proteins CheckM2 predicted on the
deposited wCchiA contigs; each query was binned by the identity and coverage of its best
hit (>=99 % over >=90 % of the query / 90-99 % / <90 % / no hit).

Run from the repository root:  python3 workflow/figures/build_figI_S6.py
Writes figures/si/FigI_S6_mag_quality.{pdf,png}.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TAB = os.path.join(REPO, "data", "tables")
OUT = os.path.join(REPO, "figures", "si")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

GREEN = "#e7f0e0"
CEIL = "#7fae6a"
BLUE, BLUE_L, BLUE_LL, GREY = "#1f6fb4", "#7fb0d6", "#c9dcec", "#d5d5d5"
ORANGE = "#d2691e"

# ---------------------------------------------------------------- input
quality = {r["Name"]: r for r in csv.DictReader(open(os.path.join(TAB, "figI2_genome_quality.csv")))}
bins = list(csv.DictReader(open(os.path.join(TAB, "figI2_removed_protein_identity.csv"))))

fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 4.3))

# ---------------------------------------------------------------- panel a
axa.add_patch(Rectangle((-0.35, 90), 5 + 0.35, 11, facecolor=GREEN, edgecolor="none", zorder=0))
axa.axvline(5, color=CEIL, lw=1.2, ls="--", zorder=1)
axa.text(4.88, 94.42, "MIMAG contamination ceiling (5 %)", color=CEIL, ha="right", va="bottom", fontsize=8.5)

# references
for x, y, lab, dx, dy, ha, va in [(0.36, 99.48, "$w$Mel (ref.)", 0.30, -0.05, "left", "center"),
                                  (0.67, 99.98, "$w$Pip (ref.)", 0.00, 0.45, "center", "bottom")]:
    axa.plot(x, y, marker="^", ms=6.5, color="#9a9a9a", ls="none", zorder=3)
    axa.text(x + dx, y + dy, lab, fontsize=8.5, ha=ha, va=va, color="#5a5a5a")

# step points of wCchiA + wCchiB
v3 = (float(quality["strainA_full"]["Contamination"]), float(quality["strainA_full"]["Completeness"]))
v4 = (float(quality["strainA_main"]["Contamination"]), float(quality["strainA_main"]["Completeness"]))
v5 = (float(quality["strainA_v5_deposited"]["Contamination"]), float(quality["strainA_v5_deposited"]["Completeness"]))
b = (float(quality["strainB"]["Contamination"]), float(quality["strainB"]["Completeness"]))

axa.plot(*v3, marker="o", ms=7, color=BLUE_LL, ls="none", zorder=3)
axa.plot(*v4, marker="o", ms=7, color=BLUE_L, ls="none", zorder=3)
axa.plot(*v5, marker="o", ms=8, color=BLUE, ls="none", zorder=4)
axa.plot(*b, marker="o", ms=8, color=ORANGE, ls="none", zorder=4)

axa.annotate("", xy=v4, xytext=v3, arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.1))
axa.annotate("", xy=v5, xytext=v4, arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.1))
axa.text(5.45, 95.55, "$w$CchiA v3\n(superseded)", fontsize=8.5, ha="left", va="center", color="#5a5a5a")
axa.text(3.05, 94.88, "remove 2 unplaced sequences\n(49,614 bp)", fontsize=8, ha="center", va="center", color="#444444")
axa.text(2.95, 95.70, "remove 2 duplicate copies\n(41,737 bp)", fontsize=8, ha="center", va="center", color="#444444")
axa.text(v5[0], v5[1] - 0.22, "$w$CchiA\n(deposited)", fontsize=8.5, ha="center", va="top", color="#14507f")
axa.text(b[0] + 0.22, b[1], "$w$CchiB", fontsize=9, ha="left", va="center", color="#8b4513")

axa.set_xlim(-0.35, 7.6)
axa.set_ylim(94.3, 101.0)
axa.set_xticks(range(0, 8))
axa.set_xlabel("CheckM2 contamination (%)")
axa.set_ylabel("CheckM2 completeness (%)")
axa.spines[["top", "right"]].set_visible(False)
axa.text(-0.13, 1.03, "a", transform=axa.transAxes, fontsize=13, fontweight="bold", va="top")

# ---------------------------------------------------------------- panel b
labels = [r["bin"].replace(">=", "$\\geq$").replace("no hit", "no hit") for r in bins]
counts = [int(r["n"]) for r in bins]
colors = [BLUE, BLUE_L, BLUE_LL, GREY]
left = 0
for lab, n, c in zip(labels, counts, colors):
    axb.barh(0, n, left=left, height=0.42, color=c, edgecolor="white", lw=0.8)
    axb.text(left + n / 2, 0, str(n), ha="center", va="center", fontsize=10,
             color="white" if c in (BLUE, BLUE_L) else "#333333")
    left += n
axb.set_xlim(0, sum(counts) + 6)
axb.set_ylim(-0.55, 0.75)
axb.set_yticks([])
axb.set_xlabel("predicted proteins on the removed sequences (n = %d)" % sum(counts))
axb.set_title("best-hit amino-acid identity to the deposited $w$CchiA assembly",
              fontsize=9.5, pad=12)
for i, (lab, n) in enumerate(zip(labels, counts)):
    centre = sum(counts[:i]) + n / 2
    axb.text(centre, -0.33, lab, ha="center", va="top", fontsize=8.5)
axb.annotate("8 proteins lack a counterpart in\nthe deposit (contig_08 3, ptg000003l 4,\nNODE_12 1); 3,538 bp of contig_08\ncarries four conserved genes",
             xy=(sum(counts) - counts[-1] / 2, 0.21), xytext=(sum(counts) - counts[-1] / 2, 0.50),
             ha="center", va="bottom", fontsize=8, color="#666666",
             arrowprops=dict(arrowstyle="-", color="#999999", lw=0.9))
axb.spines[["top", "right", "left"]].set_visible(False)
axb.tick_params(axis="y", length=0)
axb.text(-0.06, 1.02, "b", transform=axb.transAxes, fontsize=13, fontweight="bold", va="top")

fig.tight_layout()
os.makedirs(OUT, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, "FigI_S6_mag_quality." + ext), dpi=400, bbox_inches="tight")
print("wrote", os.path.join(OUT, "FigI_S6_mag_quality.{pdf,png}"))
