
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint
C_RICE, C_WO, C_COT, C_ICH, C_OTH = '#2C5F8D', '#C8792E', '#2C5F8D', '#8C8C8C', '#C8792E'
C_UNRES, C_NONCOT = '#9AA7B4', '#7FB0D8'
def wil(k, n):
    lo, hi = proportion_confint(k, n, method='wilson'); return k/n*100, lo*100, hi*100
IT = lambda s: r"$\mathit{" + s + "}$"
PAL = {"Cotesia (Braconidae)": C_COT, "Braconidae, ID unresolved": C_UNRES,
       "Braconidae, non-Cotesia": C_NONCOT, "Ichneumonidae": C_ICH, "Eupelmidae s.l.": C_OTH,
       "Braconidae, non-Cotesia/unresolved": C_NONCOT}

def build_main(mets, comp, cgroups, strat, out, dnote=None):
    fig = plt.figure(figsize=(7.09, 5.9))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.30, left=0.085, right=0.985, top=0.94, bottom=0.115)
    axA = fig.add_subplot(gs[0, 0])
    for i, (lab_, rk, rn, wk, wn, pv) in enumerate(mets):
        tops = []
        for j, (k, n, col) in enumerate([(rk, rn, C_RICE), (wk, wn, C_WO)]):
            p, lo, hi = wil(k, n); tops.append(hi)
            axA.bar(i + (j - 0.5) * 0.34, p, 0.31, color=col, edgecolor='none', zorder=2)
            axA.errorbar(i + (j - 0.5) * 0.34, p, yerr=[[p - lo], [hi - p]], fmt='none',
                         ecolor='#333333', elinewidth=0.8, capsize=1.8, zorder=3)
            axA.text(i + (j - 0.5) * 0.34, hi + 2, f"{k}/{n}", ha='center', va='bottom',
                     fontsize=5.8, color='#333333', zorder=3)
        axA.text(i, max(tops) + 9, f"$P$ = {pv:.2g}", ha='center', va='bottom', fontsize=6.2)
    axA.set_xticks(range(len(mets))); axA.set_xticklabels([m[0] for m in mets])
    axA.set_ylabel("percentage of larvae (%)"); axA.set_ylim(0, 137); axA.set_yticks([0, 25, 50, 75, 100])
    axA.legend(handles=[plt.Rectangle((0,0),1,1, color=C_RICE), plt.Rectangle((0,0),1,1, color=C_WO)],
               labels=["rice", "water oat"], frameon=False, loc='upper left', bbox_to_anchor=(-0.02, 1.02),
               handlelength=0.9, handleheight=0.9, fontsize=6.4, ncol=2, columnspacing=1.0)

    axB = fig.add_subplot(gs[0, 1])
    order = list(comp["rice"].keys())
    for i, hp in enumerate(["rice", "wateroat"]):
        tot = sum(comp[hp].values()); bot = 0
        for k in order:
            v = comp[hp][k]
            if v == 0: continue
            frac = v / tot * 100
            axB.bar(i, frac, 0.55, bottom=bot, color=PAL[k], edgecolor='white', linewidth=0.7, zorder=2)
            nm = IT("Cotesia") + "\n(Braconidae)" if k.startswith("Cotesia") else k
            if frac >= 18:
                axB.text(i, bot + frac / 2, f"{nm}\n{v}", ha='center', va='center', fontsize=6.0,
                         color='white', zorder=3, linespacing=1.35)
            else:
                side = -1 if i == 0 else 1
                axB.annotate(f"{k}  {v}", xy=(i + 0.28 * side, bot + frac / 2),
                             xytext=(i + 0.42 * side, bot + frac / 2 + (7 if i == 0 else 4)),
                             fontsize=5.6, color='#333333', va='center',
                             ha='right' if i == 0 else 'left',
                             arrowprops=dict(arrowstyle='-', lw=0.6, color='#666666'), zorder=3)
            bot += frac
        axB.text(i, 101, f"n = {tot}", ha='center', va='bottom', fontsize=6.2)
    axB.set_xticks([0, 1]); axB.set_xticklabels(["rice", "water oat"])
    axB.set_ylabel("parasitized larvae (%)"); axB.set_ylim(0, 118); axB.set_xlim(-1.15, 1.95)
    axB.set_yticks([0, 25, 50, 75, 100])

    axC = fig.add_subplot(gs[1, 0])
    rng = np.random.default_rng(0)
    xs_all = {}
    for nm, xpos, d, col, mk, fill in cgroups:
        xs_all.setdefault(xpos, []).append(d)
        axC.scatter(xpos + rng.uniform(-0.13, 0.13, len(d)), d.para_frac, s=15, marker=mk,
                    facecolor=col if fill else 'none', edgecolor=col if fill else col,
                    linewidth=0.9 if not fill else 0.35, zorder=3, label=nm if nm else None)
    for xpos, ds in xs_all.items():
        allv = np.concatenate([d.para_frac.values for d in ds])
        axC.hlines(np.median(allv), xpos - 0.28, xpos + 0.28, color='#333333', linewidth=1.1, zorder=4)
    axC.axhline(0.4, color='#999999', linestyle=(0, (3, 2)), linewidth=0.8, zorder=1)
    nx = len(set(g[1] for g in cgroups))
    axC.set_yscale('log'); axC.set_xticks(range(nx))
    axC.set_ylabel("parasitoid fraction of library"); axC.set_xlim(-0.55, nx - 0.45); axC.set_ylim(1.5e-3, 5)
    axC.legend(frameon=False, fontsize=5.3, loc='lower left', bbox_to_anchor=(-0.02, -0.02),
               handletextpad=0.2, borderpad=0.1, labelspacing=0.22, markerscale=0.9)

    axD = fig.add_subplot(gs[1, 1])
    for i, (nm, k, n, col) in enumerate(strat):
        p, lo, hi = wil(k, n)
        axD.bar(i, p, 0.5, color=col, edgecolor='none', zorder=2)
        axD.errorbar(i, p, yerr=[[p - lo], [hi - p]], fmt='none', ecolor='#333333',
                     elinewidth=0.8, capsize=2, zorder=3)
        axD.text(i, hi + 3, f"{k}/{n}", ha='center', va='bottom', fontsize=6.0, color='#333333')
    axD.set_xticks(range(len(strat))); axD.set_xticklabels([s[0] for s in strat], fontsize=5.5)
    axD.set_ylabel(IT("Wolbachia") + "+ (%)"); axD.set_ylim(0, 132); axD.set_yticks([0, 25, 50, 75, 100])
    axD.set_xlim(-0.6, len(strat) - 0.4)
    if dnote: axD.text(-0.5, 128, dnote, ha='left', va='top', fontsize=5.9)

    for ax_, L_ in [(axA, 'A'), (axB, 'B'), (axC, 'C'), (axD, 'D')]:
        ax_.text(-0.16, 1.06, L_, transform=ax_.transAxes, fontsize=10, fontweight='bold', va='bottom')
        assert not ax_.get_title()
    fig.savefig(out + ".png", dpi=400); fig.savefig(out + ".pdf")
    return fig
