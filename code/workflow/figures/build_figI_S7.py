import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

apply_figure_style(frame='open', font='Liberation Sans', sizes=(8, 7, 6))

DRY_E, DRY_F = '#2C5F8D', '#DCE7F2'
WET_E, WET_F = '#A85A1E', '#F7E6D3'
MAT_E, MAT_F = '#4A4A4A', '#EAEAEA'
SUP_E, SUP_F = '#2C5F8D', '#C8DAEC'
NEU_E, NEU_F = '#7A7A7A', '#F2F2F2'

L0, R0, W = 3.0, 53.0, 44.0
LC, RC = L0 + W / 2, R0 + W / 2

SP = r"$\mathit{C.}$ $\mathit{suppressalis}$"
GREY_T = '#444444'


def letter(x, y, s):
    ax.text(x, y, s, ha='left', va='bottom', fontsize=10, fontweight='bold', zorder=4)


def arrow(x, y_from, y_to, color='#555555'):
    ax.add_patch(FancyArrowPatch((x, y_from), (x, y_to), arrowstyle='-|>', mutation_scale=8,
                                 linewidth=0.9, color=color, shrinkA=0, shrinkB=0, zorder=1))


LH_T, LH_D, GAP_TD, PAD = 1.68, 1.52, 1.5, 1.9


def measure(title, detail):
    nT, nD = len(title.split("\n")), (len(detail.split("\n")) if detail else 0)
    c = nT * LH_T + (GAP_TD + nD * LH_D if nD else 0)
    return c + 2 * PAD, nT


def box3(x0, y1, w, title, detail, ec, fc):
    h, nT = measure(title, detail)
    ax.add_patch(FancyBboxPatch((x0, y1 - h), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                linewidth=0.9, edgecolor=ec, facecolor=fc,
                                mutation_aspect=0.5, zorder=2))
    cx, top = x0 + w / 2, y1 - PAD
    ax.text(cx, top, title, ha='center', va='top', fontsize=7, zorder=3, linespacing=1.4)
    if detail:
        ax.text(cx, top - nT * LH_T - GAP_TD, detail, ha='center', va='top', fontsize=6.2,
                color=GREY_T, zorder=3, linespacing=1.5)
    return h



SP2 = SP
dryF = [
 ("Field cohort of individual overwintering larvae",
  "247 libraries = 210 laboratory arm (paired rice / water oat)\n+ 37 public arm (external reference, PRJNA1011001)"),
 ("Read QC (fastp) and host read removal\n(" + SP + " genome)", "n = 247 carried forward"),
 ("Per-strain " + r"$\mathit{Wolbachia}$" + " mapping · parasitoid\nmitochondrial panel · kraken2 sidecar",
  "highest-breadth strain retained per library\npanel contains no " + r"$\mathit{Wolbachia}$" + " sequence"),
 ("Positivity call\nRPM \u2265 200 AND breadth \u2265 0.05; uniformity $U$ at the limit",
  r"$\mathit{Wolbachia}$+ 33/210 · parasitoid+ 45/210" "\n"
  r"0 counter-examples (Fisher $P=8.3\times10^{-29}$)" "\n"
  "parasitoid load \u2265 0.4 \u2192 29/29 " r"$\mathit{Wolbachia}$+"),
 ("Parasitoid identification\n(panel assignment + de novo mitogenomes)",
  r"$\mathit{Cotesia}$ 32/32 $\mathit{Wolbachia}$+ · Ichneumonidae 0/11" "\n"
  r"$\mathit{Cotesia}$ vs all others: Fisher $P=4.5\times10^{-10}$"),
 ("Host plant \u2192 parasitoid taxon \u2192 " + r"$\mathit{Wolbachia}$",
  r"rice 27/110 vs water oat 6/100 ($P=2.3\times10^{-4}$)" "\n"
  "parasitism rates equal ($P=0.18$); effect abolished\n"
  r"within $\mathit{Cotesia}$ (27/27 vs 5/5, $P=1.0$)"),
 ("PacBio HiFi assembly of one parasitized larva",
  "two complete genomes: " r"$w$CchiA (A) · $w$CchiB (B)" "\n"
  "ANI 85.29 % · MLST 5/5 each · MIMAG high quality\n"
  "no recombination (0/451 core genes)"),
]
wetF = [
 (r"$\mathit{Wolbachia}$-free " + SP + " colony",
  "\u2248110 generations on artificial diet\n" r"+ $\mathit{Cotesia\ chilonis}$ colony"),
 ("Experimental parasitism of the colony",
  r"larvae exposed to $\mathit{C.\ chilonis}$ females" "\ndissection confirms parasitism before assay"),
 (r"$\mathit{wsp}$ PCR (81F/691R, ≈630 bp)",
  "unparasitized larvae 0/5\nparasitized larvae 5/5\n" r"adult $\mathit{C.\ chilonis}$ 5/5"),
 ("Bidirectional Sanger sequencing\n(two adult wasps)",
  "590 / 589 bp consensus, 0 mismatches in overlap\n"
  r"100.00 % identical to the $w$CchiA $\mathit{wsp}$ allele" "\n"
  r"79.1 / 79.0 % to $w$CchiB"),
 ("FISH (Cy3\u201316S rRNA probe)",
  "signal in the wasp larva inside the host and in\nthe adult wasp abdomen; unparasitized larvae blank\n"
  "n = 3 individuals per group · no ovary/egg series"),
]
hA = [measure(t, d)[0] for t, d in dryF]
hB = [measure(t, d)[0] for t, d in wetF]
gA = (93.8 - sum(hA)) / (len(hA) - 1)
gB = min(8.0, (93.8 - sum(hB)) / (len(hB) - 1))

fig = plt.figure(figsize=(7.2, 8.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
letter(L0, 96.6, 'A'); letter(R0, 96.6, 'B')
ax.text(L0 + W/2, 96.9, 'Computational arm (247 libraries)', ha='center', va='bottom', fontsize=7, color=DRY_E)
ax.text(R0 + W/2, 96.9, 'Experimental arm (laboratory colonies)', ha='center', va='bottom', fontsize=7, color=WET_E)
y = 95.8
for i, (t, d) in enumerate(dryF):
    h = box3(L0, y, W, t, d, DRY_E, DRY_F)
    if i < len(dryF) - 1: arrow(LC, y - h, y - h - gA)
    y -= h + gA
y = 95.8
for i, (t, d) in enumerate(wetF):
    h = box3(R0, y, W, t, d, WET_E, WET_F)
    if i < len(wetF) - 1: arrow(RC, y - h, y - h - gB)
    y -= h + gB
fig.savefig("newfigs/FigI_S7_workflow.png", dpi=400)
fig.savefig("newfigs/FigI_S7_workflow.pdf")
