suppressMessages({library(readr); library(dplyr); library(ggplot2)})
set.seed(42)

d <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)
dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)

# ---- Evidence Line A: co-occurrence association ----

# Guard: Fisher.test requires a 2x2 contingency table.
# If either flag has < 2 distinct levels, skip and record NA with a note.
fisher_p   <- NA_real_
fisher_or  <- NA_real_
fisher_note <- NA_character_
ct <- table(as.integer(d$para_pos), as.integer(d$wolb_pos))
if (nrow(ct) >= 2 && ncol(ct) >= 2) {
  ft <- fisher.test(ct)
  fisher_p  <- ft$p.value
  fisher_or <- unname(ft$estimate)
} else {
  fisher_note <- paste0(
    "Fisher not computable: contingency table is ", nrow(ct), "x", ncol(ct),
    " (need 2x2). para_pos levels: {", paste(rownames(ct), collapse = ","),
    "}; wolb_pos levels: {", paste(colnames(ct), collapse = ","), "}"
  )
}

# Guard: Spearman requires non-constant vectors.
spearman_rho  <- NA_real_
spearman_p    <- NA_real_
spearman_note <- NA_character_
if (length(unique(d$wolb_rpm)) < 2 || length(unique(d$para_frac)) < 2) {
  spearman_note <- paste0(
    "Spearman not computable: zero-variance vector (",
    if (length(unique(d$wolb_rpm)) < 2)  "wolb_rpm constant",
    if (length(unique(d$wolb_rpm)) < 2 && length(unique(d$para_frac)) < 2) "; ",
    if (length(unique(d$para_frac)) < 2) "para_frac constant",
    ")"
  )
} else {
  sp <- tryCatch(
    suppressWarnings(cor.test(d$wolb_rpm, d$para_frac, method = "spearman")),
    error = function(e) NULL
  )
  if (!is.null(sp)) {
    spearman_rho <- unname(sp$estimate)
    spearman_p   <- sp$p.value
  } else {
    spearman_note <- "Spearman not computable: cor.test error (likely zero-variance or n=1)"
  }
}

stats <- tibble(
  test       = c("fisher_para_vs_wolb", "spearman_rpm_vs_parafrac"),
  statistic  = c(NA_real_, spearman_rho),
  p_value    = c(fisher_p, spearman_p),
  odds_ratio = c(fisher_or, NA_real_),
  note       = c(fisher_note, spearman_note)
)
write_tsv(stats, "results/03_profile/cooccurrence_stats.tsv")

# Figure 1: Wolbachia abundance by parasitoid status
ggsave("reports/figs/cooccurrence_box.png",
  ggplot(d, aes(factor(para_pos), log10(wolb_rpm + 1))) +
    geom_boxplot() + geom_jitter(width = .15, alpha = .5) +
    labs(x = "parasitoid positive", y = "log10(Wolbachia RPM+1)"),
  width = 5, height = 4, dpi = 150)

# ---- Evidence Line C: Wolbachia frequency stratified by host_plant ----
freq <- d %>%
  group_by(host_plant) %>%
  summarise(n = n(), wolb_rate = mean(wolb_pos), .groups = "drop")
write_tsv(freq, "results/03_profile/wolb_freq_by_host.tsv")

ggsave("reports/figs/cooccurrence_freq_host.png",
  ggplot(freq, aes(host_plant, wolb_rate)) +
    geom_col() + ylim(0, 1) +
    labs(x = "host plant", y = "Wolbachia positive rate"),
  width = 5, height = 4, dpi = 150)

# ---- Report ----
fisher_str  <- if (!is.na(fisher_p))  sprintf("p=%.3g, OR=%.3g", fisher_p, fisher_or) else fisher_note
spearman_str <- if (!is.na(spearman_rho)) sprintf("rho=%.3g, p=%.3g", spearman_rho, spearman_p) else spearman_note

report_lines <- c(
  "# Stage 3 共现分析",
  "",
  sprintf("N = %d 个样本", nrow(d)),
  "",
  "## 证据线 A：共现关联",
  "",
  sprintf("- Fisher (寄生蜂×Wolbachia): %s", fisher_str),
  sprintf("- Spearman (RPM vs 寄生蜂比例): %s", spearman_str),
  "",
  "## 证据线 C：按寄主植物分层 Wolbachia 频率",
  "",
  paste(capture.output(print(as.data.frame(freq))), collapse = "\n"),
  "",
  "见 reports/figs/cooccurrence_*.png"
)
if (is.na(fisher_p) || is.na(spearman_rho)) {
  report_lines <- c(report_lines,
    "",
    "> 注：上述检验中某一变量在当前数据中为单一类别（退化），统计量不可计算；",
    "> 这通常发生在 Wolbachia 或寄生蜂阳性样本数为 0、或全部同类时。"
  )
}
writeLines(report_lines, "reports/stage03_cooccurrence.md")

cat("done\n")
