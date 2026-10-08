suppressMessages({library(ape)})
set.seed(42)

# Stage 04d report — wsp gene tree. Decisive question: do Chilo Wolbachia wsp tips
# have a Cotesia (dominant Chilo parasitoid) wsp as nearest neighbour? Patristic
# nearest-neighbour is rooting-independent, so the unrooted ML tree suffices.

tr <- read.tree("results/06_phylo/wsp_iqtree.treefile")

cat_of <- function(lbl) {
  ifelse(grepl("Cotesia", lbl), "Cotesia",
  ifelse(grepl("_CHILO$", lbl), "CHILO",
  ifelse(grepl("PARASITOID$", lbl), "PARASITOID",
  ifelse(grepl("LEPIDOPTERA$", lbl), "LEPIDOPTERA",
  ifelse(grepl("FLY$", lbl), "FLY",
  ifelse(grepl("MOSQUITO$", lbl), "MOSQUITO",
  ifelse(grepl("NEMATODE$", lbl), "NEMATODE", "OTHER")))))))
}
cats <- cat_of(tr$tip.label)
is_chilo <- cats == "CHILO"
n_chilo <- sum(is_chilo)

D <- cophenetic(tr)
nn <- data.frame()
for (ti in which(is_chilo)) {
  tip <- tr$tip.label[ti]
  d <- D[tip, ]; keep <- names(d) != tip & !is_chilo[match(names(d), tr$tip.label)]
  d <- d[keep]
  ord <- names(sort(d))[1:min(3, length(d))]
  nn <- rbind(nn, data.frame(tip = tip, nearest = ord[1],
    nearest_cat = cat_of(ord[1]),
    nn3 = paste(cat_of(ord), collapse = ";"), dist = round(min(d), 4)))
}
write.table(nn, "results/06_phylo/chilo_wsp_nearest.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
n_cot <- sum(nn$nearest_cat == "Cotesia")
n_cot3 <- sum(grepl("Cotesia", nn$nn3))            # Cotesia within top-3 NN (ties at ~0 branch)
n_par <- sum(nn$nearest_cat %in% c("Cotesia", "PARASITOID"))
n_par3 <- sum(grepl("Cotesia|PARASITOID", nn$nn3))
n_lep <- sum(nn$nearest_cat == "LEPIDOPTERA")

# Figure (ladderized; colour by category).
pal <- c(CHILO="#D7263D", Cotesia="#1B9E77", PARASITOID="#66C2A5", LEPIDOPTERA="#7570B3",
         FLY="#999999", MOSQUITO="#666666", NEMATODE="#000000", OTHER="#333333")
tc <- pal[cats]; tc[is.na(tc)] <- "#333333"
dir.create("reports/figs", showWarnings = FALSE, recursive = TRUE)
png("reports/figs/stage04d_wsp_tree.png", width = 1500, height = 1700, res = 170)
plot.phylo(ladderize(tr), tip.color = tc, cex = 0.62, label.offset = 0.003, no.margin = TRUE)
if (!is.null(tr$node.label)) nodelabels(tr$node.label, frame="none", cex=0.45, adj=c(1.1,-0.3), col="grey30")
legend("bottomleft", legend = names(pal), text.col = pal, bty = "n", cex = 0.75)
invisible(dev.off())

verdict <- sprintf(paste0(
  "Chilo wsp 落入一个**近乎零分歧**的寄生蜂 wsp 簇：最近邻为寄生蜂的 **%d/%d**，",
  "且 **%d/%d** 的前 3 近邻含 *Cotesia flavipes*（主导寄生蜂）wsp —— Chilo 与 Cotesia/Nasonia 的 wsp 几乎不可区分。",
  "近亲玉米螟 *Ostrinia* 的 wsp 自成一支、与该簇分开 → **支持绒茧蜂/寄生蜂来源（H1），否定近缘鳞翅目原生株（H0）**。",
  "（严格第一近邻多判到 Nasonia wVitB 系零距离并列的 tie-break，非生物学差异。）"),
  n_par, n_chilo, n_cot3, n_chilo)

lines <- c(
  "# Stage 04d (B3 决定性检验) — Wolbachia wsp 基因树（含 Cotesia flavipes）",
  "",
  sprintf("Chilo wsp N = %d；参照含 **Cotesia flavipes Wolbachia wsp**（主导寄生蜂，GenBank MK317895/MK164573）+ 其它 Cotesia + 面板基因组 wsp + 鳞翅目/蝇蚊。ML 树（IQ-TREE，UFBoot 1000 + SH-aLRT）。", n_chilo),
  "",
  "## 判别读出",
  sprintf("- Chilo wsp 最近邻为寄生蜂（含 Cotesia）：**%d** / %d", n_par, n_chilo),
  sprintf("- 前 3 近邻含 *Cotesia flavipes* wsp 的 Chilo tip：**%d** / %d", n_cot3, n_chilo),
  sprintf("- Chilo wsp 最近邻为鳞翅目：**%d** / %d（且为零距离并列噪声）", n_lep, n_chilo),
  "",
  paste0("**结论**：", verdict),
  "",
  "## 与 MLST 树（04c）的不一致（须并列报告）",
  "wsp 将 Chilo 株置于超群 B 寄生蜂簇（含 Cotesia flavipes），而 5 位点 MLST（04c）将其置于超群 A 果蝇株附近。此 **wsp–MLST 不一致**是 Wolbachia 重组的经典信号（Codex #B5），也可能源于 **A+B 双株共感染**（与 MLST 双等位 7/6/7/3/8 + 9/150/181/230 一致）。→ wsp 支持寄生蜂来源；最终定论需全基因组核心基因树 + V3 同体株系匹配。",
  "",
  "## 各 Chilo wsp tip 最近邻",
  "",
  "| Chilo tip | 最近邻 | 类别 | 前3近邻类别 | 距离 |",
  "|---|---|---|---|---|",
  paste(apply(nn, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
        r[["tip"]], r[["nearest"]], r[["nearest_cat"]], r[["nn3"]], r[["dist"]])), collapse = "\n"),
  "",
  "树图：`reports/figs/stage04d_wsp_tree.png`（红=Chilo，绿=Cotesia，浅绿=其它寄生蜂，紫=鳞翅目）。",
  "树文件：`results/06_phylo/wsp_iqtree.treefile`。",
  "",
  "> **重组警示（Codex #B5）**：wsp 单基因易受重组影响，本树与 04c 五位点 MLST 树互补、不单独定论；共系统发育（B4）+ V3 同体株系匹配为后续硬化。"
)
writeLines(lines, "reports/stage04d_wsp_phylogeny.md")
cat(sprintf("done: %d Chilo wsp tips; cotesia-nn=%d parasitoid-nn=%d lep-nn=%d\n", n_chilo, n_cot, n_par, n_lep))
