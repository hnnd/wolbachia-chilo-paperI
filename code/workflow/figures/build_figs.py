"""Rebuild the 01.wolb manuscript figures with pubfig (publication-chart-skill workflow).

Data: figdata/results/** pulled from hpc-c. Panels are composed by passing matplotlib
axes into pubfig plot calls, then exported through pubfig's Nature figure spec.
Phylogenetic trees (old Fig4A/B) have no pubfig chart family and are NOT rebuilt here.
"""
import dataclasses
import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pubfig as pf
from scipy import stats

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# --- theme: keep the Nature spec but use an installed metric-compatible face ---
_base = pf.get_theme("nature")
THEME = dataclasses.replace(_base, font_family=["Liberation Sans", "DejaVu Sans", "sans-serif"])
pf.register_theme("nature-local", THEME)
pf.set_default_theme("nature-local")

D = "figdata/results"
prof = pd.read_csv(f"{D}/03_profile/per_sample_profile.tsv", sep="\t")
carr = pd.read_csv(f"{D}/06_denovo/fig4_carriage.tsv", sep="\t", header=None,
                   names=["sample", "family", "para_reads", "wolb_rpm", "wolb_pos"])
ab = pd.read_csv(f"{D}/07_cophylo/wolb_ab.tsv", sep="\t")
coef = pd.read_csv(f"{D}/03_profile/multivariate_coef.tsv", sep="\t")
sweep = pd.read_csv(f"{D}/03_profile/threshold_sweep.tsv", sep="\t")

assert len(prof) == 247 and int(prof.wolb_pos.sum()) == 31 and int(prof.para_pos.sum()) == 45


def letter(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")


def rate_bar(ax, rates, labels, ylab, counts):
    """Proportion bars with n/N annotated above each bar.

    Label x-positions are read from the drawn bar geometry, NOT from the
    category index: pubfig lays categories out at i * category_spacing
    (0.75 by default), so integer positions drift right and the last label
    can fall outside the axes.
    """
    pf.bar(np.asarray(rates), category_names=list(labels), x_label="",
           y_label=ylab, legend_show=False, ax=ax)
    ax.set_ylim(0, 1.12 if max(rates) > 0.5 else max(rates) * 1.35)
    bars = [p for p in ax.patches]
    assert len(bars) == len(rates), f"expected {len(rates)} bars, drew {len(bars)}"
    pad = ax.get_ylim()[1] * 0.02
    for p, r, n in zip(bars, rates, counts):
        xc = p.get_x() + p.get_width() / 2
        assert abs(p.get_height() - r) < 1e-9, "bar height/rate mismatch"
        ax.text(xc, r + pad, n, ha="center", va="bottom", fontsize=6)


out = {}

# ============================ Figure 1 — Evidence A ============================
fig1, ax1 = plt.subplots(1, 3, figsize=(7.2, 2.4))

ct = pd.crosstab(prof.wolb_pos.astype(int), prof.para_pos.astype(int))
mat = np.array([[ct.loc[1, 1], ct.loc[1, 0]], [ct.loc[0, 1], ct.loc[0, 0]]], float)
pf.heatmap(mat, category_names=["Parasitoid +", "Parasitoid \u2212"], annotate=True,
           annotate_fmt=".0f", colorscale="Blues", cbar=False,
           imshow_aspect_square="auto", ax=ax1[0])
ax1[0].set_yticks([0, 1]); ax1[0].set_yticklabels(["Wolbachia +", "Wolbachia \u2212"])
# tick labels already carry both factors; drop the confusion-matrix defaults
ax1[0].set_xlabel(""); ax1[0].set_ylabel("")
# no titles anywhere on the figures — cohort size / subset go in the caption
letter(ax1[0], "A")

pf.scatter(prof.para_frac.values, prof.wolb_breadth.values,
           labels=np.where(prof.para_pos == 1, "Parasitoid +", "Parasitoid \u2212"),
           x_label="Parasitoid read fraction", y_label="Wolbachia genome breadth",
           ax=ax1[1])
ax1[1].axhline(0.05, color="0.55", ls=":", lw=0.6)
letter(ax1[1], "B")

hp = prof.groupby("host_plant").agg(n=("wolb_pos", "size"), pos=("wolb_pos", "sum"))
hp = hp.reindex(["rice", "wateroat"])
rate_bar(ax1[2], (hp.pos / hp.n).values, ["Rice", "Water-oat"],
         "Wolbachia+ frequency", [f"{int(p)}/{int(n)}" for p, n in zip(hp.pos, hp.n)])
letter(ax1[2], "C")

fig1.tight_layout()
out["Fig1_cooccurrence"] = fig1

# ============================ Figure 2 — dose & sufficiency ====================
fig2, ax2 = plt.subplots(1, 2, figsize=(6.0, 2.5))

sub = prof[prof.para_pos == 1]
rho, pv = stats.spearmanr(sub.para_frac, sub.wolb_breadth)
pf.scatter(sub.para_frac.values, sub.wolb_breadth.values,
           x_label="Parasitoid read fraction", y_label="Wolbachia genome breadth",
           show_regression=True, legend_show=False, ax=ax2[0])
letter(ax2[0], "A")

bins = [-0.001, 0.001, 0.01, 0.1, 0.4, 1.01]
labs = ["<0.001", "0.001\u20130.01", "0.01\u20130.1", "0.1\u20130.4", "\u22650.4"]
b = prof.assign(bin=pd.cut(prof.para_frac, bins=bins, labels=labs)) \
        .groupby("bin", observed=True).agg(n=("wolb_pos", "size"), pos=("wolb_pos", "sum"))
rate_bar(ax2[1], (b.pos / b.n).values, labs,
         "Wolbachia+ frequency", [f"{int(p)}/{int(n)}" for p, n in zip(b.pos, b.n)])
ax2[1].set_xlabel("Parasitoid read fraction")
ax2[1].tick_params(axis="x", labelrotation=30)
letter(ax2[1], "B")

fig2.tight_layout()
out["Fig2_dose"] = fig2

# ============================ Figure 3 — Evidence D ============================
fig3, ax3 = plt.subplots(1, 2, figsize=(6.4, 2.5))

pf.scatter(np.log10(carr.para_reads.values),
           np.log10(carr.wolb_rpm.values + 0.1),
           labels=carr.family.values,
           x_label="log$_{10}$ parasitoid reads",
           y_label="log$_{10}$ Wolbachia RPM", ax=ax3[0])
letter(ax3[0], "A")

fam_order = ["Braconidae (Cotesia)", "Tachinidae (fly)", "Ichneumonidae"]
f = carr.groupby("family").agg(n=("wolb_pos", "size"), pos=("wolb_pos", "sum")).reindex(fam_order)
rate_bar(ax3[1], (f.pos / f.n).values, ["Braconidae\n(Cotesia)", "Tachinidae", "Ichneumonidae"],
         "Wolbachia+ frequency", [f"{int(p)}/{int(n)}" for p, n in zip(f.pos, f.n)])
letter(ax3[1], "B")

fig3.tight_layout()
out["Fig3_taxon_carriage"] = fig3

# ============================ Figure 4 — Evidence B (non-tree) =================
fig4, ax4 = plt.subplots(1, 2, figsize=(6.2, 2.5))

pf.scatter(ab.depthA.values, ab.depthB.values,
           x_label="Supergroup A depth (\u00d7)", y_label="Supergroup B depth (\u00d7)",
           legend_show=False, ax=ax4[0])
_lim = [0, max(ab.depthA.max(), ab.depthB.max()) * 1.08]
ax4[0].plot(_lim, _lim, ls="--", lw=0.6, color="0.65", zorder=0)
ax4[0].text(0.96, 0.90, "y = x", transform=ax4[0].transAxes, ha="right",
            va="top", fontsize=6, color="0.5")
letter(ax4[0], "A")

order = ab.sort_values("Bfrac")
one = [pf.get_palette("nature")[0]] * len(order)
pf.bar(order.Bfrac.values, category_names=list(order["sample"]), x_label="",
       y_label="Supergroup B fraction", legend_show=False,
       color_palette=one, ax=ax4[1])
ax4[1].axhline(ab.Bfrac.mean(), color="0.25", ls="--", lw=0.7)
ax4[1].set_ylim(0, 0.62)
ax4[1].text(0.02, ab.Bfrac.mean() + 0.015, f"mean {ab.Bfrac.mean():.2f}",
            transform=ax4[1].get_yaxis_transform(), ha="left", va="bottom",
            fontsize=6, color="0.25")
ax4[1].tick_params(axis="x", labelrotation=90)
letter(ax4[1], "B")

fig4.tight_layout()
out["Fig4_evidenceB_panels"] = fig4

# ============================ Figure S2 — Firth ORs ============================
c = coef[coef.term != "(Intercept)"].copy()
name_map = {"z_lpf": "Parasitoid load (z)",
            "z_depth": "Sequencing depth (z)",
            "z_biomass": "Non-host biomass (z)"}
c["label"] = c.term.map(name_map)
figS2 = pf.forest_plot(c.OR.values, c.lo.values, c.hi.values, labels=list(c.label),
                       x_label="Odds ratio (95 % CI)", reference=1.0, x_scale="log")
out["FigS2_multivariate_or"] = figS2

# ============================ Figure S3 — threshold sweep ======================
piv = sweep.pivot_table(index="b_thr", columns="p_thr", values="counterex", aggfunc="max")
figS3 = pf.heatmap(piv.values.astype(float), annotate=True, annotate_fmt=".0f",
                   colorscale="Reds", zmin=0, zmax=2,
                   cbar_label="Wolbachia+ / parasitoid\u2212 counter-examples",
                   x_label="Parasitoid fraction threshold",
                   y_label="Wolbachia breadth threshold")
axS3 = figS3.axes[0]
axS3.set_xticks(range(len(piv.columns))); axS3.set_xticklabels([f"{v:g}" for v in piv.columns])
axS3.set_yticks(range(len(piv.index))); axS3.set_yticklabels([f"{v:g}" for v in piv.index])
out["FigS3_threshold_sensitivity"] = figS3

# ==================== Figure S1 — cohort overview =============================
# (A) sampling design: samples per collection site, split by host plant
# (B) sequencing depth by cohort arm
# (C) marker positivity by cohort arm — the public arm is deeper yet wholly negative
figS1 = plt.figure(figsize=(13.2, 4.2))
gsS1 = figS1.add_gridspec(1, 3, width_ratios=[1.5, 1.0, 1.0], wspace=0.42)
axA, axB, axC = (figS1.add_subplot(gsS1[0, i]) for i in range(3))

# One public library (SAMN01797446) has no recorded site; groupby would drop it
# silently and the panel would show 246 of 247 samples. Give it an explicit bin.
_geo = prof.geo.fillna("Public, site n.r.")
_site = (prof.assign(geo=_geo).groupby(["geo", "host_plant"]).size().unstack(fill_value=0)
         .reindex(columns=["rice", "wateroat"], fill_value=0))
_site = _site.loc[_site.sum(axis=1).sort_values().index]
assert _site.values.sum() == len(prof), f"panel A plots {_site.values.sum()} of {len(prof)}"
pf.bar(_site.values.astype(float), category_names=list(_site.index),
       series_names=["Rice", "Water-oat"], orientation="horizontal",
       x_label="", y_label="Samples", legend_show=True, ax=axA)
axA.set_ylabel("")          # site codes are self-labelling; avoids a clipped axis title
letter(axA, "A")

_depth = [(prof.loc[prof.source == s, "total"] / 1e6).values for s in ["lab", "sra"]]
pf.box(_depth, category_names=["Field-collected", "Public data"], ax=axB)
axB.set_ylabel("Sequencing depth (M read pairs)")
axB.set_yscale("log")
letter(axB, "B")

_arm = prof.groupby("source").agg(n=("wolb_pos", "size"), wolb=("wolb_pos", "sum"),
                                  para=("para_pos", "sum")).loc[["lab", "sra"]]
_ar = np.column_stack([(_arm.wolb / _arm.n).values, (_arm.para / _arm.n).values])
pf.bar(_ar, category_names=["Field-collected", "Public data"],
       series_names=["Wolbachia +", "Parasitoid +"],
       x_label="", y_label="Positivity rate", ax=axC)
axC.set_ylim(0, 0.30)
_exp1 = [_ar[0, 0], _ar[0, 1], _ar[1, 0], _ar[1, 1]]
_cnt1 = [f"{int(_arm.wolb.iloc[0])}/{int(_arm.n.iloc[0])}", f"{int(_arm.para.iloc[0])}/{int(_arm.n.iloc[0])}",
         f"{int(_arm.wolb.iloc[1])}/{int(_arm.n.iloc[1])}", f"{int(_arm.para.iloc[1])}/{int(_arm.n.iloc[1])}"]
_b1 = sorted(axC.patches, key=lambda q: q.get_x())
assert len(_b1) == 4, f"expected 4 grouped bars, drew {len(_b1)}"
for _q, _v, _c in zip(_b1, _exp1, _cnt1):
    assert abs(_q.get_height() - _v) < 1e-9, "S1C bar order mismatch"
    axC.text(_q.get_x() + _q.get_width() / 2, _v + 0.006, _c, ha="center", va="bottom", fontsize=6)
letter(axC, "C")
out["FigS1_cohort_overview"] = figS1

# ==================== Figure S4 — host-plant stratification ===================
# Both rates side by side: parasitism prevalence is essentially flat across hosts
# while Wolbachia carriage is not — the difference tracks parasitoid TAXON
# (water-oat is Ichneumonidae-dominated, rice is Cotesia-dominated), not parasitism rate.
# LAB-COLLECTED ONLY. The 37 public (SRA) datasets are all rice and all negative for
# both markers, so pooling them into the rice denominator dilutes rice parasitism from
# 25.5% to 19.0% and manufactures an apparent "equal parasitism" between hosts.
# Only the lab cohort samples both hosts under one regime, so only it is comparable.
_lab = prof[prof.source == "lab"]
assert len(_lab) == 210 and prof[prof.source == "sra"].para_pos.sum() == 0
_hp = _lab.groupby("host_plant").agg(n=("wolb_pos", "size"), wolb=("wolb_pos", "sum"),
                                     para=("para_pos", "sum")).loc[["rice", "wateroat"]]
# pubfig reads data[i][j] as category i, series j (NOT series x category)
_rates = np.column_stack([(_hp.wolb / _hp.n).values, (_hp.para / _hp.n).values])
figS4 = pf.bar(_rates, category_names=["Rice", "Water-oat"],
               series_names=["Wolbachia +", "Parasitoid +"],
               x_label="Host plant", y_label="Positivity rate")
axS4 = figS4.axes[0]
axS4.set_ylim(0, 0.26)
_expected = [_rates[0, 0], _rates[0, 1], _rates[1, 0], _rates[1, 1]]
_counts = [f"{int(_hp.wolb.iloc[0])}/{int(_hp.n.iloc[0])}", f"{int(_hp.para.iloc[0])}/{int(_hp.n.iloc[0])}",
           f"{int(_hp.wolb.iloc[1])}/{int(_hp.n.iloc[1])}", f"{int(_hp.para.iloc[1])}/{int(_hp.n.iloc[1])}"]
_bars = sorted(axS4.patches, key=lambda p: p.get_x())
assert len(_bars) == 4, f"expected 4 grouped bars, drew {len(_bars)}"
for _p, _v, _c in zip(_bars, _expected, _counts):
    assert abs(_p.get_height() - _v) < 1e-9, f"grouped bar order mismatch: {_p.get_height()} vs {_v}"
    axS4.text(_p.get_x() + _p.get_width() / 2, _v + 0.005, _c,
              ha="center", va="bottom", fontsize=6)
out["FigS4_host_plant"] = figS4

# ---- guard: no titles are ever drawn on the figures (they belong in captions) ----
for _name, _fig in out.items():
    _titles = [a.get_title() for a in _fig.axes if a.get_title()]
    _sup = _fig._suptitle.get_text() if _fig._suptitle else ""
    assert not _titles and not _sup, f"{_name} has in-figure title(s): {_titles} {_sup!r}"

# ============================ export ==========================================
paths = {}
for name, fig in out.items():
    w = "single" if name in ("FigS3_threshold_sensitivity", "FigS4_host_plant") else "double"
    paths[name] = pf.batch_export(fig, f"newfigs/{name}", formats=("png", "pdf"),
                                  spec="nature", width=w, dpi=300)

STATS = {
    "n": len(prof), "wolb_pos": int(prof.wolb_pos.sum()), "para_pos": int(prof.para_pos.sum()),
    "contingency_wolbpos_parapos": int(ct.loc[1, 1]), "counterexamples": int(ct.loc[1, 0]),
    "spearman_rho": round(float(rho), 3), "spearman_p": float(pv),
    "rice_rate": round(float(hp.loc["rice", "pos"] / hp.loc["rice", "n"]), 4),
    "wateroat_rate": round(float(hp.loc["wateroat", "pos"] / hp.loc["wateroat", "n"]), 4),
    "highload_bin": f"{int(b.pos.iloc[-1])}/{int(b.n.iloc[-1])}",
    "carriage": {k: f"{int(v.pos)}/{int(v.n)}" for k, v in f.iterrows()},
    "ichneumonidae_max_reads": int(carr[carr.family == "Ichneumonidae"].para_reads.max()),
    "Bfrac_mean": round(float(ab.Bfrac.mean()), 3), "Bfrac_sd": round(float(ab.Bfrac.std()), 3),
    "sweep_rows": len(sweep),
    "sweep_zero_counterex_at_b_ge_0.05": int((sweep[sweep.b_thr >= 0.05].counterex == 0).all()),
    "sweep_nonzero_rows": int((sweep.counterex > 0).sum()),
}
