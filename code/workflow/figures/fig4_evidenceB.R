#!/usr/bin/env Rscript
# Fig 4 — Evidence line B: Wolbachia strain phylogeny + A+B co-infection.
# No titles are drawn on the figure (figure/panel titles belong in the caption).
# Panels: (A) wsp ML tree, (B) MLST ML tree, (C) dual-supergroup co-infection
# scatter, (D) B4 A:B dosage invariance. Sources: Stage 04c/04d/04f/04g.
suppressMessages({ library(ggtree); library(ggplot2); library(patchwork); library(ape) })
have_midpoint <- requireNamespace("phangorn", quietly = TRUE)

figdir <- "reports/figs"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---- Okabe-Ito colorblind-safe palette by lineage (FIXED order) -----------
LINEAGES <- c("Chilo (this study)","Cotesia (parasitoid)","Other parasitoid",
              "Lepidoptera","Drosophila","Mosquito","Nematode")
pal <- setNames(c("#D55E00","#0072B2","#56B4E9","#009E73","#E69F00","#CC79A7","#555555"), LINEAGES)
classify <- function(lab) {
  g <- ifelse(grepl("_CHILO$|^Chilo", lab), "Chilo (this study)",
       ifelse(grepl("Cotesia", lab), "Cotesia (parasitoid)",
       ifelse(grepl("_PARASITOID|wLcla|wVitB|wEnc|wTpre|wTkay|wUni", lab), "Other parasitoid",
       ifelse(grepl("_LEPIDOPTERA|wEel|wOfur|wOscap", lab), "Lepidoptera",
       ifelse(grepl("_FLY|wNo|wRi|wMel", lab), "Drosophila",
       ifelse(grepl("_MOSQUITO|wPip", lab), "Mosquito", "Nematode"))))))
  factor(g, levels = LINEAGES)
}
# reference strain -> clean display name (host in parentheses)
REFMAP <- c(
  wCflav_MK164573.1="wCflav (C. flavipes)", "wCflav_CF_MK317895.1"="wCflav-CF (C. flavipes)",
  wCmel_OR597565.1="wCmel (Cotesia)", wCgl_PQ472555.1="wCgl (C. glomerata)",
  wLcla="wLcla (Leptopilina)", wVitB="wVitB (Nasonia)", wEnc="wEnc (Encarsia)",
  wTpre="wTpre (Trichogramma)", wTkay="wTkay (Trichogramma)", wUni="wUni (Muscidifurax)",
  wEel="wEel (Ephestia)", wOfur="wOfur (Ostrinia)", wOscap="wOscap (Ostrinia)",
  wNo="wNo (Drosophila)", wRi="wRi (Drosophila)", wMel="wMel (Drosophila)",
  wPip="wPip (Culex)", wBm="wBm (Brugia, outgroup)")
prettytip <- function(lab) {
  out <- character(length(lab))
  for (i in seq_along(lab)) {
    l <- lab[i]
    if (grepl("_CHILO$", l)) { out[i] <- gsub("_CHILO$", "", gsub("^Chilo_", "", l)); next }
    key <- l
    key <- gsub("^[0-9]+_PARASITOID_", "", key)          # drop leading "1_PARASITOID_"
    key <- gsub("_Cotesia$", "", key)                     # drop trailing "_Cotesia"
    key <- gsub("_(PARASITOID|LEPIDOPTERA|FLY|MOSQUITO|NEMATODE)$", "", key)
    out[i] <- if (!is.na(REFMAP[key])) REFMAP[key] else gsub("_", " ", key)
  }
  out
}

tree_panel <- function(nwk, xaxis_lab, midpoint = FALSE) {
  tr <- read.tree(nwk)
  if (midpoint && have_midpoint) tr <- phangorn::midpoint(tr)
  tr <- ladderize(tr)
  grp <- classify(tr$tip.label)
  dd <- data.frame(label = tr$tip.label, lineage = grp, disp = prettytip(tr$tip.label))
  p <- ggtree(tr, size = 0.4) %<+% dd +
    geom_tippoint(aes(color = lineage), size = 1.8) +
    geom_tiplab(aes(label = disp, color = lineage), size = 2.1, offset = 0.002, show.legend = FALSE) +
    scale_color_manual(values = pal, drop = FALSE, name = "Lineage", limits = LINEAGES) +
    theme_tree2() +
    xlab(xaxis_lab) +
    theme(axis.title.x = element_text(size = 9),
          legend.position = "bottom", legend.title = element_text(size=9),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, label = "")))
  p + xlim(0, max(p$data$x, na.rm = TRUE) * 1.35)
}

# ---- A: wsp tree (decisive) ; B: MLST tree --------------------------------
pA <- tree_panel("results/06_phylo/wsp_iqtree.treefile",
                 "wsp tree \u2014 substitutions per site (midpoint-rooted)",
                 midpoint = TRUE)
pB <- tree_panel("results/06_phylo/iqtree.treefile",
                 "MLST (5-locus) tree \u2014 substitutions per site") +
  guides(color = "none")   # single shared lineage legend comes from panel A

# ---- C: co-infection dual-supergroup scatter (04f) -------------------------
co <- read.delim("results/07_cophylo/fig3_coinfect.tsv", header = FALSE,
                 col.names = c("sample","breadthA","breadthB","het","xgenome","call"))
co$het <- as.numeric(co$het); co$xgenome <- as.numeric(co$xgenome)
pC <- ggplot(co, aes(breadthA, breadthB)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey70") +
  geom_point(aes(size = xgenome, fill = het*100), shape = 21, color = "grey20", alpha = .9) +
  scale_fill_viridis_c(option = "C", name = "het %") +
  scale_size_continuous(range = c(2.5, 8), name = "assembly\n(x genome)") +
  labs(x = "breadth on supergroup A (wMel)",
       y = "breadth on supergroup B (wPip)") +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0.5,1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.99,0.01), legend.justification = c(1,0),
        legend.background = element_rect(fill=alpha("white",.7), color=NA),
        legend.key.size = unit(0.35,"cm"))

# ---- D: B4 A:B dosage invariance (04g) ------------------------------------
ab <- read.delim("results/07_cophylo/wolb_ab.tsv")
ab$Bfrac <- as.numeric(ab$Bfrac)
mB <- mean(ab$Bfrac); sdB <- sd(ab$Bfrac)
ab <- ab[order(ab$Bfrac), ]; ab$sample <- factor(ab$sample, levels = ab$sample)
pD <- ggplot(ab, aes(sample, Bfrac)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=mB-sdB, ymax=mB+sdB, fill="#0072B2", alpha=.12) +
  geom_hline(yintercept = mB, linetype = 2, color = "#0072B2") +
  geom_point(size = 3, shape = 21, fill = "#0072B2", color = "grey20") +
  annotate("text", x = 1, y = 0.62, hjust = 0, size = 3,
           label = sprintf("B-fraction = %.3f ± %.3f (SD)\ninvariant across samples", mB, sdB)) +
  labs(x = NULL, y = "Wolbachia B-fraction  depthB/(A+B)") +
  coord_cartesian(ylim = c(0.25, 0.65)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &   # no figure title: it belongs in the caption
  theme(plot.tag = element_text(face = "bold", size = 13), legend.position = "bottom")

ggsave(file.path(figdir, "Fig4_evidenceB.png"), fig, width = 12, height = 11, dpi = 300)
ggsave(file.path(figdir, "Fig4_evidenceB.pdf"), fig, width = 12, height = 11)
cat("[Fig4] wrote", file.path(figdir, "Fig4_evidenceB.png/.pdf"), "\n")
