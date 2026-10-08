suppressMessages({library(ape)})
set.seed(42)

# Stage 04c (B3) report — read the MLST ML tree and quantify whether the Chilo
# Wolbachia tips nest among PARASITOID Wolbachia (H1) or with Lepidoptera Wolbachia
# (H0). Host type is the last underscore field of each tip label.

tr <- read.tree("results/06_phylo/iqtree.treefile")
htype <- sub(".*_", "", tr$tip.label)            # CHILO / PARASITOID / LEPIDOPTERA / FLY / MOSQUITO / NEMATODE
is_chilo <- htype == "CHILO"
n_chilo <- sum(is_chilo)

# Nearest-neighbour host type for each Chilo tip (patristic distance, k=3 non-Chilo).
D <- cophenetic(tr)
nn_tab <- data.frame()
for (ti in which(is_chilo)) {
  tip <- tr$tip.label[ti]
  d <- D[tip, ]
  d <- d[names(d) != tip & !is_chilo[match(names(d), tr$tip.label)]]
  ord <- names(sort(d))[1:min(3, length(d))]
  nnh <- sub(".*_", "", ord)
  nn_tab <- rbind(nn_tab, data.frame(
    tip = tip,
    nearest = ord[1], nearest_host = nnh[1],
    nn3 = paste(nnh, collapse = ";"),
    nn_dist = round(min(d), 4)))
}
maj <- names(sort(table(nn_tab$nearest_host), decreasing = TRUE))[1]
n_para <- sum(nn_tab$nearest_host == "PARASITOID")
n_lep  <- sum(nn_tab$nearest_host == "LEPIDOPTERA")
write.table(nn_tab, "results/06_phylo/chilo_nearest_neighbour.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figure: tree coloured by host type, Chilo tips bold/red.
pal <- c(CHILO = "#D7263D", PARASITOID = "#1B9E77", LEPIDOPTERA = "#7570B3",
         FLY = "#999999", MOSQUITO = "#666666", NEMATODE = "#000000")
tip_col <- pal[htype]; tip_col[is.na(tip_col)] <- "#333333"
dir.create("reports/figs", showWarnings = FALSE, recursive = TRUE)
png("reports/figs/stage04c_mlst_tree.png", width = 1500, height = 1700, res = 170)
plot.phylo(tr, tip.color = tip_col, cex = 0.7, label.offset = 0.002, no.margin = TRUE)
if (!is.null(tr$node.label))
  nodelabels(tr$node.label, frame = "none", cex = 0.45, adj = c(1.1, -0.3), col = "grey30")
legend("bottomleft", legend = names(pal), text.col = pal, bty = "n", cex = 0.8)
invisible(dev.off())

verdict <- paste0(
  "5 位点 MLST 树**不解析/倾向超群 A**：多数 Chilo tip 最近邻为果蝇株（wRi/wMel，与 best-strain 比对一致），",
  "**0 个**最近邻为寄生蜂。**但此结论受限**：(i) 主导寄生蜂 *Cotesia* 的 Wolbachia 无基因组、未入本面板（缺最关键参照）；",
  "(ii) MLST 仅 5 位点、骨架支持率低、多数样本仅部分位点回收；(iii) 近亲玉米螟 *Ostrinia* 亦**不**与 Chilo 聚类。",
  "→ MLST 本身不足以判别，须看 **wsp 树（Stage 04d，含 Cotesia flavipes）**：wsp 将 Chilo 置于寄生蜂簇。",
  "两标记不一致 = Wolbachia 重组 / 共感染信号（Codex #B5）。")

lines <- c(
  "# Stage 04c (B3) — Wolbachia MLST 株系系统发育",
  "",
  sprintf("纳入 Chilo 样本 Wolbachia 序列 N = %d（≥3 个 MLST 位点回收），与寄主标注参照面板（寄生蜂 / 鳞翅目 / 蝇蚊 / 外群）联合建 5 位点分区 ML 树（IQ-TREE，UFBoot 1000 + SH-aLRT，外群 wBm/Brugia）。", n_chilo),
  "",
  "## 判别读出",
  "",
  sprintf("- Chilo tip 最近邻为寄生蜂：**%d** / %d", n_para, n_chilo),
  sprintf("- Chilo tip 最近邻为鳞翅目：**%d** / %d", n_lep, n_chilo),
  sprintf("- 多数最近邻 host_type：**%s**", maj),
  "",
  paste0("**结论**：", verdict),
  "",
  "## 各 Chilo tip 最近邻（patristic）",
  "",
  "| Chilo tip | 最近邻 | host | 前3近邻 host | 距离 |",
  "|---|---|---|---|---|",
  paste(apply(nn_tab, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
        r[["tip"]], r[["nearest"]], r[["nearest_host"]], r[["nn3"]], r[["nn_dist"]])), collapse = "\n"),
  "",
  "树图：`reports/figs/stage04c_mlst_tree.png`（红=Chilo，绿=寄生蜂，紫=鳞翅目）。",
  "树文件：`results/06_phylo/iqtree.treefile`；近邻表：`results/06_phylo/chilo_nearest_neighbour.tsv`。",
  "",
  "> 注意：MLST 易受重组影响（Codex #B5）；本树为 5 位点串联分区 ML，wsp 与全基因组核心基因树及共系统发育（B4）为后续加固。低覆盖样本仅部分位点回收，缺失位点以 gap 补齐。"
)
writeLines(lines, "reports/stage04c_phylogeny.md")
cat(sprintf("done: %d Chilo tips, parasitoid-nn=%d lepidoptera-nn=%d\n", n_chilo, n_para, n_lep))
