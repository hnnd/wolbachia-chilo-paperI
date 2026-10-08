import numpy as np, pandas as pd, dataclasses, logging
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pubfig as pf
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
_base = pf.get_theme("nature")
THEME = dataclasses.replace(_base, font_family=["Liberation Sans","DejaVu Sans","sans-serif"])
pf.register_theme("nature-local", THEME); pf.set_default_theme("nature-local")
prof = pd.read_csv("hpc/scp-90a348a1b44a872d/pd/per_sample_profile.tsv", sep="\t")
_u = pd.read_csv("work/coverage_uniformity_all247.csv")[["sample","wolb_pos_rev"]]
prof = prof.merge(_u, on="sample", how="left")
prof["wolb_pos"] = np.where(prof.wolb_pos_rev.notna(), prof.wolb_pos_rev, prof.wolb_pos)
assert len(prof) == 247 and int(prof.wolb_pos.sum()) == 33 and int(prof.para_pos.sum()) == 45


def letter(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")


import os
os.makedirs("newfigs", exist_ok=True)

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

# ---- guard: no titles are ever drawn on the figures (they belong in captions) ----
_titles = [a.get_title() for a in figS1.axes if a.get_title()]
_sup = figS1._suptitle.get_text() if figS1._suptitle else ""
assert not _titles and not _sup, f"FigI_S1_cohort_overview has in-figure title(s): {_titles} {_sup!r}"

paths = pf.batch_export(figS1, "newfigs/FigI_S1_cohort_overview", formats=("png", "pdf"),
                        spec="nature", width="double", dpi=300)

