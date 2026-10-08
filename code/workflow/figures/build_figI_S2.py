import os
os.makedirs('newfigs', exist_ok=True)
import tarfile, glob, os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

apply_figure_style(frame='open', font='Liberation Sans', sizes=(8, 7, 6.4))

# Load the diagnostic data from the two tarballs
os.makedirs("diag", exist_ok=True)
with tarfile.open("hpc/scp-90a348a1b44a872d/diagtar/wolb_diag.tgz") as tf:
    tf.extractall("diag")

NMC = ["sample", "ref", "len", "nm", "ident", "pos", "mapq"]
nm = pd.concat([pd.read_csv(f, sep="\t", names=NMC) for f in sorted(glob.glob("diag/*.nm.tsv"))
                if not os.path.basename(f).startswith("SAMN")])
cov = pd.concat([pd.read_csv(f, sep="\t", names=["sample", "ref", "pos", "depth"])
                 for f in sorted(glob.glob("diag/*.cov.tsv"))
                 if not os.path.basename(f).startswith("SAMN")])
nm = nm[nm.ref == "NC_002978.6"]
cov = cov[cov.ref == "NC_002978.6"]

ORD = ["NL_R_6", "NC2025_R_3", "PYH_R_3", "AQ_R_2", "HZ_R_10", "NL_W_14"]

G_WMEL = 1267782

def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum())

rows = []
for s in ORD:
    d = cov[cov["sample"] == s]
    if not len(d):
        continue
    tot = d.depth.sum()
    dec = d.nlargest(max(1, len(d) // 10), "depth").depth.sum()
    rows.append(dict(sample=s, cov_sites=len(d), mean_depth_covered=round(d.depth.mean(), 3),
                     max_depth=int(d.depth.max()), top_decile_share=round(dec / tot, 3),
                     gini=round(gini(d.depth.values), 3),
                     span_kb=round((d.pos.max() - d.pos.min()) / 1000, 1),
                     n_windows_10kb=d.pos.floordiv(10000).nunique(),
                     pct_genome_windows=round(d.pos.floordiv(10000).nunique() / (G_WMEL // 10000) * 100, 1)))
cs = pd.DataFrame(rows)

# Load the public arm diagnostic data
with tarfile.open("hpc/scp-90a348a1b44a872d/diagtar/wolb_diag2.tgz") as tf:
    tf.extractall("diag")

nm2 = pd.concat([pd.read_csv(f, sep="\t", names=NMC) for f in sorted(glob.glob("diag/SAMN*.nm.tsv"))])
cov2 = pd.concat([pd.read_csv(f, sep="\t", names=["sample", "ref", "pos", "depth"])
                  for f in sorted(glob.glob("diag/SAMN*.cov.tsv"))])
nm2w = nm2[nm2.ref == "NC_002978.6"]
cov2w = cov2[cov2.ref == "NC_002978.6"]
PUB = sorted(cov2w["sample"].unique())

pub_rows = []
for s in PUB:
    d = cov2w[cov2w["sample"] == s]
    n = len(nm2w[nm2w["sample"] == s])
    tot = d.depth.sum()
    dec = d.nlargest(max(1, len(d) // 10), "depth").depth.sum()
    pub_rows.append(dict(sample=s, cls="public over-RPM", reads=n, cov_sites=len(d),
                         sites_per_read=round(len(d) / n, 2), max_depth=int(d.depth.max()),
                         top_decile_share=round(dec / tot, 3), gini=round(gini(d.depth.values), 3),
                         n_windows_10kb=d.pos.floordiv(10000).nunique(),
                         pct_genome_windows=round(d.pos.floordiv(10000).nunique() / (G_WMEL // 10000) * 100, 1),
                         ident_median=round(nm2w[nm2w["sample"] == s].ident.median(), 4)))

cs["sites_per_read"] = [round(cs.cov_sites.iloc[i] / nm[nm["sample"] == cs["sample"].iloc[i]].shape[0], 2)
                        for i in range(len(cs))]
cs["ident_median"] = [round(nm[nm["sample"] == s].ident.median(), 4) for s in cs["sample"]]
cs["cls"] = ["reclassified positive", "reclassified positive", "noise call",
             "confirmed positive", "confirmed positive", "confirmed positive"]
cs["reads"] = [nm[nm["sample"] == s].shape[0] for s in cs["sample"]]

allcal = pd.concat([cs, pd.DataFrame(pub_rows)], ignore_index=True)

nmA = pd.concat([nm, nm2w])
covA = pd.concat([cov, cov2w])

CLS = {"NL_R_6": "r", "NC2025_R_3": "r", "AQ_R_2": "p", "HZ_R_10": "p", "NL_W_14": "p",
       "PYH_R_3": "n", **{s: "u" for s in PUB}}
CC = {"r": '#B03030', "p": '#2C5F8D', "n": '#8C8C8C', "u": '#7B4FA0'}
CLSNAME = {"r": "reclassified positive (third criterion)", "p": "confirmed positive (both criteria)",
           "n": "below both criteria", "u": "public arm, over-RPM"}
SO = ["NL_R_6", "NC2025_R_3", "AQ_R_2", "HZ_R_10", "NL_W_14", "PYH_R_3"] + PUB
spr = allcal.set_index("sample")

fig = plt.figure(figsize=(7.09, 8.2))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.35, 0.95], hspace=0.36, wspace=0.42,
                      left=0.185, right=0.985, top=0.975, bottom=0.065)

axA = fig.add_subplot(gs[0, 0])
data = [nmA[nmA["sample"] == s].ident.values for s in SO]
bp = axA.boxplot(data, vert=False, widths=0.62, patch_artist=True, showfliers=False,
                 medianprops=dict(color='#222222', linewidth=0.9),
                 whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7))
for patch, s in zip(bp['boxes'], SO):
    patch.set(facecolor=CC[CLS[s]], edgecolor='none', alpha=0.85)
axA.set_yticks(range(1, len(SO) + 1))
axA.set_yticklabels(SO, fontsize=5.7)
axA.invert_yaxis()
axA.set_xlim(0.80, 1.01)
axA.set_xlabel("per-read alignment identity to wMel")

axB = fig.add_subplot(gs[0, 1])
for i, s in enumerate(SO):
    v = spr.loc[s, "sites_per_read"]
    axB.barh(i, v, 0.62, color=CC[CLS[s]], edgecolor='none', zorder=2)
    axB.text(v * 1.25, i, f"{v:.1f}", va='center', fontsize=5.7, color='#333333')
axB.set_xscale('log')
axB.set_xlim(0.4, 700)
axB.set_yticks(range(len(SO)))
axB.set_yticklabels([])
axB.invert_yaxis()
axB.set_xlabel("covered sites per aligned read")
axB.axvline(50, color='#666666', linestyle=(0, (3, 2)), linewidth=0.8, zorder=1)
axB.text(52, len(SO) - 0.4, "50", fontsize=5.7, color='#666666', va='bottom')

axC = fig.add_subplot(gs[1, :])
for i, s in enumerate(SO):
    d = covA[covA["sample"] == s]
    axC.scatter(d.pos.values / 1e6, np.full(len(d), i), s=0.8, marker='|',
                color=CC[CLS[s]], alpha=0.5, linewidth=0.35, zorder=3)
axC.set_yticks(range(len(SO)))
axC.set_yticklabels([f"{s}   ({int(spr.loc[s, 'reads']):,} reads, "
                     f"{spr.loc[s, 'pct_genome_windows']:.0f} % windows)" for s in SO], fontsize=5.7)
axC.invert_yaxis()
axC.set_xlim(0, G_WMEL / 1e6)
axC.set_ylim(len(SO) - 0.4, -0.6)
axC.set_xlabel("position on the wMel genome (Mb)")
axC.legend(handles=[plt.Line2D([], [], color=CC[k], lw=3) for k in ["p", "r", "n", "u"]],
           labels=[CLSNAME[k] for k in ["p", "r", "n", "u"]], frameon=False, fontsize=5.6,
           loc='lower left', bbox_to_anchor=(0.005, -0.015), handlelength=1.2, labelspacing=0.25,
           ncol=2, columnspacing=1.2)


# ---- panel D: depth profile, pile-up artefact versus genuine infection -------
axD = fig.add_subplot(gs[2, :])
WIN = 2000
for s, lab, col in [("SAMN10449037", "SAMN10449037 (public arm, over-RPM)", CC["u"]),
                    ("AQ_R_2", "AQ_R_2 (confirmed positive)", CC["p"])]:
    d = covA[covA["sample"] == s]
    w = (d.pos // WIN).astype(int)
    prof = d.groupby(w).depth.mean()
    x = prof.index.values * WIN / 1e6
    axD.vlines(x, 0.55, prof.values, color=col, lw=0.55, alpha=0.85,
               label=f"{lab}: max {int(d.depth.max()):,}x, top-decile share {spr.loc[s,'top_decile_share']:.2f}")
axD.set_yscale("log")
axD.set_ylim(0.55, 3000)
axD.set_xlim(0, G_WMEL / 1e6)
axD.set_xlabel("position on the wMel genome (Mb)")
axD.set_ylabel("mean depth in\n2-kb window (x)")
axD.legend(frameon=False, fontsize=5.6, loc="upper left", handlelength=1.2, labelspacing=0.25)

for ax_, L_ in [(axA, 'A'), (axB, 'B'), (axC, 'C'), (axD, 'D')]:
    ax_.text({axA:-0.30, axB:-0.06}.get(ax_, -0.255), 1.02,
             L_, transform=ax_.transAxes, fontsize=10, fontweight='bold', va='bottom')
    assert not ax_.get_title()

fig.savefig("newfigs/FigI_S2_signal_authenticity.png", dpi=400)
fig.savefig("newfigs/FigI_S2_signal_authenticity.pdf")