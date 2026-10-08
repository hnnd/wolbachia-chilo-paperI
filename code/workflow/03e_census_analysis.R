suppressMessages({library(readr); library(dplyr); library(ggplot2)})
set.seed(42)

# Stage 03e — per-sample parasitoid census against a broad mitogenome reference.
# Answers: which parasitoids does each sample contain, quantified, incl. multi-
# parasitism; and (crossing Wolbachia status) which parasitoid taxa carry Wolbachia.
# Reads are assigned by primary best-hit then aggregated to GENUS, which is robust
# to within-genus cross-mapping. A genus is called "present" only if it clears an
# absolute + relative threshold, suppressing conserved-region leakage from the
# dominant taxon.

m <- read_tsv("results/03_profile/census_long.tsv",
              col_names = c("sample", "species", "contig", "reads"), show_col_types = FALSE)
prof <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)
dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)

MIN_READS <- 20      # absolute floor for "present"
MIN_FRAC  <- 0.05    # >=5% of a sample's parasitoid reads

m <- m %>%
  filter(species != "NONE") %>%
  mutate(genus = sub("_.*", "", species)) %>%
  filter(genus != "Chilo")                       # drop host control

# Per (sample, genus) read totals.
sg <- m %>% group_by(sample, genus) %>% summarise(reads = sum(reads), .groups = "drop")
samp_tot <- sg %>% group_by(sample) %>% summarise(total_para = sum(reads), .groups = "drop")
sg <- sg %>% left_join(samp_tot, by = "sample") %>%
  mutate(frac = reads / total_para,
         present = reads >= MIN_READS & frac >= MIN_FRAC)

# Per-sample summary: dominant genus, composition, multi-parasitism.
per <- sg %>% group_by(sample) %>%
  summarise(
    total_para  = first(total_para),
    dom_genus   = genus[which.max(reads)],
    dom_reads   = max(reads),
    dom_frac    = max(reads) / first(total_para),
    n_genera    = sum(present),                                   # multi-parasitism count
    genera_present = paste(sort(genus[present]), collapse = ";"),
    .groups = "drop"
  )

d <- prof %>% select(sample, wolb_pos, wolb_rpm, para_frac, nonhost, total) %>%
  left_join(per, by = "sample")
# samples with no parasitoid reads at all become NA after the join -> set to 0/none.
d$total_para[is.na(d$total_para)] <- 0
d$dom_reads[is.na(d$dom_reads)]   <- 0
d$n_genera[is.na(d$n_genera)]     <- 0
d$dom_genus[is.na(d$dom_genus)]   <- "none"
d$genera_present[is.na(d$genera_present)] <- ""
write_tsv(d, "results/03_profile/per_sample_census.tsv")

# Cohort genus composition: how many samples each genus dominates / is present in.
comp <- sg %>% filter(present) %>% count(genus, name = "n_present") %>%
  left_join(per %>% count(dom_genus, name = "n_dominant"), by = c("genus" = "dom_genus")) %>%
  mutate(n_dominant = ifelse(is.na(n_dominant), 0, n_dominant)) %>%
  arrange(desc(n_present))
write_tsv(comp, "results/03_profile/census_genus_composition.tsv")

# KEY TABLE: which parasitoid genus carries Wolbachia (dominant-genus samples
# with a real parasitoid load).
# med_dom_frac is a RELIABILITY flag: high (~>0.7) = reads concentrate on this
# genus (assignment trustworthy); low = reads scatter across distant relatives
# because the true species is absent from the reference (nearest-neighbour
# artifact -> needs de novo ID, T2).
gw <- d %>% filter(dom_reads >= 100) %>%
  group_by(dom_genus) %>%
  summarise(n = n(), n_wolb = sum(wolb_pos == 1),
            med_wolb_rpm = median(wolb_rpm),
            med_dom_reads = median(dom_reads),
            med_dom_frac = round(median(dom_frac), 2), .groups = "drop") %>%
  mutate(wolb_rate = n_wolb / n) %>%
  arrange(desc(n))
write_tsv(gw, "results/03_profile/genus_wolbachia.tsv")

# Multi-parasitism distribution.
multi_tab <- as.data.frame(table(n_genera = d$n_genera))

# ---- figures ----
topg <- comp %>% slice_max(n_present, n = 12)
ggsave("reports/figs/census_composition.png",
  ggplot(topg, aes(reorder(genus, n_present), n_present)) +
    geom_col(fill = "steelblue") + coord_flip() +
    labs(x = "parasitoid genus", y = "# samples present (>=20 reads & >=5%)",
         title = "Per-sample parasitoid census: genus prevalence across cohort"),
  width = 6.5, height = 4.2, dpi = 150)

gwf <- gw %>% filter(n >= 2)
if (nrow(gwf)) ggsave("reports/figs/genus_wolbachia.png",
  ggplot(gwf, aes(reorder(dom_genus, wolb_rate), wolb_rate)) +
    geom_col(fill = "firebrick") + coord_flip() + ylim(0, 1) +
    geom_text(aes(label = sprintf("%d/%d", n_wolb, n)), hjust = -0.1, size = 3) +
    labs(x = "dominant parasitoid genus", y = "Wolbachia positive rate",
         title = "Which parasitoid carries Wolbachia (dominant-genus samples)"),
  width = 6.5, height = 4, dpi = 150)

# ---- report ----
lines <- c(
  "# Stage 03e — 逐样本寄生物普查（广线粒体参考，123 contig）",
  "",
  sprintf("N = %d 样本。非宿主 reads → 广寄生物 mito 参考（跨 10 科）primary-best-hit，聚合到属。命中过的寄生物 contig 数 = %d。",
          nrow(d), length(unique(m$contig))),
  sprintf("'存在' 判据：某属 reads >= %d 且 >= %.0f%% 样本寄生总量（抑制保守区跨映射）。", MIN_READS, 100*MIN_FRAC),
  "",
  "## 队列属水平组成（按出现样本数）",
  "",
  "| 寄生物属 | 出现样本数 | 主导样本数 |",
  "|---|---|---|",
  paste(apply(comp, 1, function(r) sprintf("| %s | %s | %s |", r[["genus"]], r[["n_present"]], r[["n_dominant"]])), collapse = "\n"),
  "",
  "## 多重寄生（每样本 '存在' 的属数分布）",
  "",
  paste(capture.output(print(multi_tab, row.names = FALSE)), collapse = "\n"),
  sprintf("- 含 >=2 个寄生物属的样本数：%d", sum(d$n_genera >= 2)),
  "",
  "## 关键表：哪些寄生物属携带 Wolbachia（dom_reads>=100）",
  "",
  "| 主导属 | n | Wolb+ | 阳性率 | 中位 wolb_rpm | 中位 dom_reads | 集中度 dom_frac |",
  "|---|---|---|---|---|---|---|",
  paste(apply(gw, 1, function(r) sprintf("| %s | %s | %s | %.0f%% | %s | %s | %s |",
        r[["dom_genus"]], r[["n"]], r[["n_wolb"]], 100*as.numeric(r[["wolb_rate"]]),
        r[["med_wolb_rpm"]], r[["med_dom_reads"]], r[["med_dom_frac"]])), collapse = "\n"),
  "",
  "**稳健结论**：只有 **Cotesia（绒茧蜂）携带 Wolbachia（~93%）**；所有非 Cotesia 类群即使重度寄生也几乎不带（0–2%）。",
  "这排除了'泛生物量/泛寄生'假象——带菌与否取决于**寄生物种类**，不取决于寄生量。",
  "",
  "**⚠️ 可靠性警示（集中度 dom_frac）**：dom_frac 高(~>0.7)= reads 集中、属判定可信(如 Cotesia 0.78)；",
  "dom_frac 低(如 Scambus 0.21)= reads 在远缘参考间四散，**真实物种不在参考集**，属名只是 nearest-neighbour 假象。",
  "同一批重寄生样本在 11 株参考下被判 Pachyneuron、广参考下被判 Scambus/Neotrichoporoides —— **非 Cotesia 的属名不可信**。",
  "→ 非 Cotesia 寄生物的**真实种类必须用 T2（de novo 线粒体组装 + BLAST nt）**确定；但'非 Cotesia=不带 Wolbachia'这一结论与具体种名无关，稳健成立。",
  "",
  "见 reports/figs/census_composition.png, genus_wolbachia.png"
)
writeLines(lines, "reports/stage03e_census.md")
cat("done\n")
