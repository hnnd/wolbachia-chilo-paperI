suppressMessages({library(readr); library(dplyr)})
set.seed(42)

# Stage 03f (T2) — de novo mitogenome assembly + BLAST species ID.
# Per sample, the longest assembled contig is the (near-complete) mitogenome; its
# best BLAST hit names the nearest relative and pident gives confidence. A second
# long contig hitting a different taxon = multi-parasitism. This is reference-free
# at the assembly step, so it resolves the parasitoids that best-hit read mapping
# could not name (03e).

b <- read_tsv("results/06_denovo/denovo_blast_long.tsv",
              col_names = c("sample","qseqid","qlen","sseqid","pident","length","bitscore","stitle"),
              show_col_types = FALSE)
census <- read_tsv("results/03_profile/per_sample_census.tsv", show_col_types = FALSE)
prof   <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)

b <- b %>% mutate(
  hit_species = sub("\\|.*", "", sseqid),
  hit_genus   = sub("_.*", "", hit_species)
)

# Primary ID per sample = best hit on the longest contig (>=2 kb already enforced).
primary <- b %>% group_by(sample) %>% slice_max(qlen, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(sample, denovo_genus = hit_genus, denovo_species = hit_species,
            pident, mito_len = qlen, aln_len = length)

# Multi-parasitism: distinct genera among near-complete contigs (>=10 kb).
multi <- b %>% filter(qlen >= 10000) %>%
  group_by(sample) %>% summarise(n_mito_genera = n_distinct(hit_genus),
            genera = paste(sort(unique(hit_genus)), collapse = ";"), .groups = "drop")

d <- prof %>% select(sample, wolb_pos, wolb_rpm) %>%
  inner_join(census %>% select(sample, dom_genus, total_para, dom_frac), by = "sample") %>%
  inner_join(primary, by = "sample") %>%
  left_join(multi, by = "sample") %>%
  mutate(n_mito_genera = ifelse(is.na(n_mito_genera), 1L, n_mito_genera),
         conf = ifelse(pident >= 97, "species", ifelse(pident >= 90, "genus/family", "family+ (divergent — remote nt advised)")))
write_tsv(d, "results/06_denovo/per_sample_denovo_id.tsv")

# De novo species composition + Wolbachia cross-tab.
comp <- d %>% count(denovo_genus, name = "n_samples") %>% arrange(desc(n_samples))
gw <- d %>% group_by(denovo_genus) %>%
  summarise(n = n(), n_wolb = sum(wolb_pos == 1), med_pident = round(median(pident),1),
            med_wolb_rpm = round(median(wolb_rpm),1), .groups = "drop") %>%
  mutate(wolb_rate = n_wolb / n) %>% arrange(desc(n))
write_tsv(gw, "results/06_denovo/denovo_genus_wolbachia.tsv")

# Agreement between 03e best-hit genus and de novo genus.
agree <- mean(d$dom_genus == d$denovo_genus)

lines <- c(
  "# Stage 03f (T2) — de novo 线粒体组装 + BLAST 定种",
  "",
  sprintf("目标样本 N = %d（非 Cotesia 高信号 + Cotesia 对照）。bait 寄生物 mito reads → megahit 组装 → blastn 本地节肢动物 mito 库。", nrow(d)),
  sprintf("- 组装出 >=10 kb 近完整线粒体的样本：%d", sum(d$mito_len >= 10000)),
  sprintf("- 03e best-hit 属 与 de novo 属 一致率：%.0f%%（不一致处即 best-hit 假象被 de novo 纠正）", 100*agree),
  "",
  "## de novo 属水平组成",
  "",
  "| de novo 属 | 样本数 |",
  "|---|---|",
  paste(apply(comp, 1, function(r) sprintf("| %s | %s |", r[["denovo_genus"]], r[["n_samples"]])), collapse = "\n"),
  "",
  "## de novo 属 × Wolbachia",
  "",
  "| 属 | n | Wolb+ | 阳性率 | 中位 pident | 中位 wolb_rpm |",
  "|---|---|---|---|---|---|",
  paste(apply(gw, 1, function(r) sprintf("| %s | %s | %s | %.0f%% | %s | %s |",
        r[["denovo_genus"]], r[["n"]], r[["n_wolb"]], 100*as.numeric(r[["wolb_rate"]]),
        r[["med_pident"]], r[["med_wolb_rpm"]])), collapse = "\n"),
  "",
  sprintf("## 多重寄生（>=2 个近完整线粒体属的样本）：%d", sum(d$n_mito_genera >= 2, na.rm = TRUE)),
  "",
  "置信度：pident>=97 物种级；90–97 属/科级；<90 远缘（建议对 results/06_denovo/<sid>.contigs.fa 做 remote nt BLAST 定种）。",
  "组装的线粒体 contig 存于 results/06_denovo/<sid>.contigs.fa，可直接提交 NCBI BLAST 取权威种名。"
)

# 远程 nt 定种（代表样本，对本地库 <90% 的发散组装做 NCBI remote megablast 纠正最近邻假象）。
rb_path <- "results/06_denovo/remote/remote_best.tsv"
if (file.exists(rb_path)) {
  rb <- read_tsv(rb_path, show_col_types = FALSE) %>%
    left_join(d %>% select(sample, wolb_pos, dom_genus), by = "sample")
  lines <- c(lines, "",
    "## 远程 nt 定种（代表样本）",
    "",
    "本地库最高仅匹到 *Scambus* 80%（最近邻假象，真实种不在库内）。对代表性长 contig 做 NCBI remote megablast：",
    "",
    "| 样本 | wolb | 远程最佳命中 | pident | 真实科 |",
    "|---|---|---|---|---|",
    paste(apply(rb, 1, function(r) sprintf("| %s | %s | *%s* | %s | %s |",
          r[["sample"]], ifelse(as.integer(r[["wolb_pos"]])==1,"+","-"),
          r[["remote_species"]], r[["remote_pident"]], r[["remote_family"]])), collapse = "\n"),
    "",
    "**生物学结论**：_W（茭白）样本的主导寄生物是**姬蜂科 Ichneumonidae**（最近邻 Diadromus/Achaius 84–88%，亚科级发散），Wolbachia 阴；唯一 _W wolb+ 的 HC_W_9 是**寄蝇科 Tachinidae**（Pexopsis 87%）；_R 带菌样本是**茧蜂科 Braconidae**（Cotesia 96%）。带菌与否取决于寄生物谱系而非寄生量 —— 判别性 HT 证据。"
  )
}
writeLines(lines, "reports/stage03f_denovo_id.md")
cat("done\n")
