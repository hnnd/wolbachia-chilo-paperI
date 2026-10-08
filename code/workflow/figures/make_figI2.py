import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from Bio import Phylo
import re

# Load data
T = Phylo.read('/data2/user/wangys/.claude-science/orgs/b0ca9587-c297-462d-8e57-4e8cf9533311/artifacts/proj_edea35f7ab19/907ed4a8-18cf-4f40-bea2-ff93c4d17271/v053b3de4_core_iqtree.contree', "newick")
T.root_with_outgroup({"name": "wBm"})

A = pd.read_csv('/data2/user/wangys/.claude-science/orgs/b0ca9587-c297-462d-8e57-4e8cf9533311/artifacts/proj_edea35f7ab19/461b335e-7012-4501-97c7-56aa8ea692f0/v16ffdab1_ani_all.tsv', sep="\t", header=None,
                names=["q","s","ani","fm","ft"])
clean = lambda x: re.sub(r".*/g_|\.fna$", "", x)
A["q"] = A.q.map(clean); A["s"] = A.s.map(clean)
_NM = {"strainA": "wCchiA", "strainB": "wCchiB"}
A["q"] = A.q.replace(_NM); A["s"] = A.s.replace(_NM)

PAF = pd.read_csv('/data2/user/wangys/.claude-science/orgs/b0ca9587-c297-462d-8e57-4e8cf9533311/artifacts/proj_edea35f7ab19/bee61eba-bb96-4509-887b-51319e03a98a/v0333edaa_AvsB.paf', sep="\t", header=None,
                  names=["q","qlen","qs","qe","strand","t","tlen","ts","te","matches","alnlen"],
                  usecols=range(11))

GCt = pd.read_csv('/data2/user/wangys/.claude-science/orgs/b0ca9587-c297-462d-8e57-4e8cf9533311/artifacts/proj_edea35f7ab19/e6891875-6c37-4daf-b62e-a5014e201315/ve9704e7a_strain_gene_content.tsv', sep="\t").set_index("metric").value

# Rename tips
NAME = {"strainA": "wCchiA", "strainB": "wCchiB"}
for t in T.get_terminals():
    if t.name in NAME: t.name = NAME[t.name]

SGA, SGB = ["wMel", "wRi", "wUni"], ["wPip", "wVitB", "wNo", "wTpre"]
SG = {**{k: "A" for k in ["wCchiA"] + SGA}, **{k: "B" for k in ["wCchiB"] + SGB}, "wBm": "D"}
ORD = ["wCchiA","wMel","wRi","wUni","wCchiB","wVitB","wPip","wNo","wTpre","wBm"]
C_A, C_B, GREY = '#2C5F8D', '#B03030', '#8C8C8C'
COL = {"A": C_A, "B": C_B, "D": GREY}

# Build ANI matrix
M = pd.DataFrame(np.nan, index=ORD, columns=ORD)
for _, rr in A.iterrows():
    if rr.q in ORD and rr.s in ORD: M.loc[rr.q, rr.s] = rr.ani
sym2 = (M + M.T) / 2
np.fill_diagonal(sym2.values, 100.0)
sym2 = sym2.combine_first(M)

isfocal = lambda n: n.startswith("wCchi")

def layout(tree):
    tips = tree.get_terminals()
    yp = {id(t): float(i) for i, t in enumerate(tips)}
    def yof(cl):
        if id(cl) in yp: return yp[id(cl)]
        v = np.mean([yof(ch) for ch in cl.clades]); yp[id(cl)] = v; return v
    yof(tree.root)
    xp = {}
    def xof(cl, acc=0.0):
        xp[id(cl)] = acc
        for ch in cl.clades: xof(ch, acc + (ch.branch_length or 0.0))
    xof(tree.root)
    return yp, xp

lanes = [("M", None), ("colony\n(no wasp)", "0/5"), ("parasitized\nlarvae", "5/5"),
         ("adult\n$\\mathit{C.\\ chilonis}$", "5/5"), ("NTC", "0")]

fig = plt.figure(figsize=(7.09, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.02], hspace=0.30, wspace=0.34,
                      left=0.075, right=0.975, top=0.965, bottom=0.075)

# Panel A: core-gene tree
axA = fig.add_subplot(gs[0, 0])
yp, xp = layout(T); xm = max(xp.values())
for cl in T.get_nonterminals():
    ys = [yp[id(ch)] for ch in cl.clades]
    axA.plot([xp[id(cl)]]*2, [min(ys), max(ys)], color='#333333', lw=0.85, solid_capstyle='round', zorder=2)
    for ch in cl.clades:
        axA.plot([xp[id(cl)], xp[id(ch)]], [yp[id(ch)]]*2, color='#333333', lw=0.85,
                 solid_capstyle='round', zorder=2)
    if cl.confidence is not None and cl is not T.root and len(cl.get_terminals()) > 1:
        axA.text(xp[id(cl)] - 0.008*xm, yp[id(cl)] - 0.30, f"{int(cl.confidence)}",
                 fontsize=5.4, color='#333333', ha='right', va='center', zorder=4)
for t in T.get_terminals():
    bold = isfocal(t.name)
    axA.text(xp[id(t)] + 0.035*xm, yp[id(t)], t.name, fontsize=6.6, va='center', ha='left',
             color=COL[SG[t.name]], fontweight=('bold' if bold else 'normal'), zorder=4)
    if bold:
        axA.scatter(xp[id(t)], yp[id(t)], s=14, marker='o', facecolor=COL[SG[t.name]],
                    edgecolor='none', zorder=5)
axA.set_ylim(len(T.get_terminals()) - 0.5, -0.9); axA.set_xlim(-0.02*xm, xm*1.62)
for sp in ("top","right","left"): axA.spines[sp].set_visible(False)
axA.set_yticks([]); axA.set_xlabel("amino-acid substitutions per site")
axA.tick_params(axis='x', labelsize=5.8)
axA.text(0.985, 0.03, "630 single-copy core genes\n226,156 aa · Q.BIRD+F+I+G4 · UFBoot 1000",
         transform=axA.transAxes, fontsize=5.5, ha='right', va='bottom', color='#555555', linespacing=1.4)

# Panel B: ANI heatmap
axB = fig.add_subplot(gs[0, 1])
Mv = sym2.loc[ORD, ORD].values.astype(float)
im = axB.imshow(Mv, cmap="viridis", vmin=80, vmax=100, aspect='equal')
axB.set_xticks(range(len(ORD))); axB.set_yticks(range(len(ORD)))
axB.set_xticklabels(ORD, rotation=90, fontsize=5.5); axB.set_yticklabels(ORD, fontsize=5.5)
for lbls in (axB.get_xticklabels(), axB.get_yticklabels()):
    for tl in lbls:
        tl.set_color(COL[SG[tl.get_text()]])
        if isfocal(tl.get_text()): tl.set_fontweight('bold')
for i in range(len(ORD)):
    for j in range(len(ORD)):
        axB.text(j, i, f"{Mv[i,j]:.1f}", ha='center', va='center', fontsize=4.1,
                 color=('white' if Mv[i,j] < 93 else 'black'))
iA, iB = ORD.index("wCchiA"), ORD.index("wCchiB")
for (rw, cl_) in [(iA, iB), (iB, iA)]:
    axB.add_patch(Rectangle((cl_-0.5, rw-0.5), 1, 1, fill=False, edgecolor='#FFFFFF', lw=1.3, zorder=5))
for sp in axB.spines.values(): sp.set_visible(False)
cb = fig.colorbar(im, ax=axB, fraction=0.042, pad=0.03)
cb.set_label("fastANI (%)", fontsize=6); cb.ax.tick_params(labelsize=5.5); cb.outline.set_visible(False)

# Panel C: synteny dotplot
axC = fig.add_subplot(gs[1, 0])
for _, rr in PAF.iterrows():
    xs = [rr.ts/1e6, rr.te/1e6]
    ys = [rr.qs/1e6, rr.qe/1e6] if rr.strand == "+" else [rr.qe/1e6, rr.qs/1e6]
    axC.plot(xs, ys, color=('#333333' if rr.strand == "+" else '#B03030'), lw=0.7,
             alpha=0.85, solid_capstyle='round')
axC.set_xlabel("wCchiA (Mb)"); axC.set_ylabel("wCchiB (Mb)")
axC.set_xlim(0, 1.40); axC.set_ylim(0, 1.40)
axC.plot([0, 1.40], [0, 1.40], color='#CCCCCC', lw=0.6, ls=(0, (4, 3)), zorder=1)
axC.legend(handles=[Line2D([], [], color='#333333', lw=1.1, label="same strand"),
                    Line2D([], [], color='#B03030', lw=1.1, label="inverted")],
           frameon=True, framealpha=0.92, edgecolor='none', fontsize=5.8, loc='upper left',
           handlelength=1.4, handletextpad=0.4, borderpad=0.3)


axD = fig.add_subplot(gs[1, 1])
_vmax = float(np.percentile(body, 99.8))
axD.imshow(body, cmap="gray", aspect="auto", vmin=0.0, vmax=_vmax, interpolation="bilinear")
axD.set_xticks([]); axD.set_yticks([])
for _s in axD.spines.values(): _s.set_visible(False)
H, Wd = body.shape
MK = [(66, "2000"), (130, "1000"), (161, "750"), (200, "500"), (266, "250"), (325, "100")]
for _y, _lab in MK:
    axD.text(-12, _y, _lab, fontsize=5.0, ha="right", va="center", color="#333333")
axD.text(-12, -20, "bp", fontsize=5.0, ha="right", va="bottom", color="#333333")
LX = {0: 55.5}
for _i in range(1, 16): LX[_i] = 110 + (_i - 1) * 54.7
axD.text(LX[0], -6, "M", fontsize=5.4, ha="center", va="bottom", color="#333333")
for _i in range(1, 16):
    axD.text(LX[_i], -6, str(_i), fontsize=5.0, ha="center", va="bottom", color="#333333")
for _a, _b, _lab in [(1, 5, "parasitised\nlarvae"), (6, 10, "adult wasps"),
                     (11, 15, "unparasitised\nlarvae")]:
    _x0, _x1 = LX[_a] - 24, LX[_b] + 24
    axD.plot([_x0, _x1], [H + 14, H + 14], color="#333333", lw=0.9, clip_on=False, zorder=4)
    axD.text((_x0 + _x1) / 2, H + 24, _lab, fontsize=5.0, ha="center", va="top",
             color="#333333", linespacing=1.3)
axD.annotate("~630 bp", xy=(LX[10] + 28, 168), xytext=(LX[13], 168), fontsize=5.2,
             color="#333333", va="center", ha="center",
             arrowprops=dict(arrowstyle="->", lw=0.7, color="#333333"))
axD.set_xlim(-38, Wd + 4)
axD.set_ylim(H + 62, -38)

for ax_, L_, dx in [(axA, 'A', -0.055), (axB, 'B', -0.20), (axC, 'C', -0.135), (axD, 'D', -0.045)]:
    ax_.text(dx, 1.01, L_, transform=ax_.transAxes, fontsize=10, fontweight='bold', va='bottom')
    assert not ax_.get_title()

fig.savefig("FigI2_strain_genomes.png", dpi=400)