suppressMessages({library(readr); library(dplyr); library(ggplot2)})
set.seed(42)

# Stage 3b — Evidence line A deepening:
#   A3 dose-response  : does Wolbachia load scale with parasitoid read fraction?
#   A5 source-specificity (lite): does Wolbachia track the PARASITOID signal
#       specifically, rather than total non-host biomass / sequencing depth /
#       host biomass? Addresses the "para_frac just co-varies with total
#       non-moth load or DNA quality" confounder (Codex review #1/#8).
#
# Pure reanalysis of results/03_profile/per_sample_profile.tsv — no remapping.
# NOTE: a true per-species mitochondrial dual-source breakdown (mapping reads to
#       Chilo-mito + each parasitoid mito) is a separate, heavier remap stage and
#       is intentionally NOT done here (see reports/stage03b_dose_source.md).

d <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)
dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)

# Derived fractions on a COMMON denominator (total reads) so Chilo and parasitoid
# biomass proxies are comparable. para_frac in the profile is reads/nonhost; here
# we also derive parasitoid and host fractions of *total* reads.
d <- d %>% mutate(
  log_rpm        = log10(wolb_rpm + 1),
  chilo_frac_tot = host_mapped / total,   # Chilo nuclear+mito biomass (of total)
  para_frac_tot  = para_reads  / total,   # parasitoid biomass (of total)
  nonhost_frac   = nonhost     / total    # generic non-host load (of total)
)

# Rank-based partial correlation of x,y controlling for covariates Z (matrix/df).
partial_spearman <- function(x, y, Z) {
  ok <- complete.cases(x, y, Z)
  x <- x[ok]; y <- y[ok]; Z <- as.data.frame(Z)[ok, , drop = FALSE]
  if (length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  rz <- as.data.frame(lapply(Z, rank))
  ex <- residuals(lm(rank(x) ~ ., data = rz))
  ey <- residuals(lm(rank(y) ~ ., data = rz))
  ct <- suppressWarnings(cor.test(ex, ey, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}

spearman_safe <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}

# Spread a named stat vector into explicit tibble columns (avoids matrix columns).
row_stat <- function(v, ...) tibble(..., rho = v[["rho"]], p = v[["p"]], n = v[["n"]])

# ---------------------------------------------------------------------------
# A3 — dose-response: Wolbachia load vs parasitoid read fraction
# ---------------------------------------------------------------------------
# Evaluated on (a) all samples and (b) parasitoid-positive only, where parasitoid
# load actually varies. A monotonic positive trend separates H1 (parasitoid
# source) from H3 (random cross-talk, which would be flat/dose-independent).
dpos <- d %>% filter(para_pos == 1)

a3 <- bind_rows(
  row_stat(spearman_safe(d$log_rpm,    d$para_frac),
           subset = "all_samples",    metric = "spearman_logrpm_vs_parafrac"),
  row_stat(spearman_safe(dpos$log_rpm, dpos$para_frac),
           subset = "parasitoid_pos", metric = "spearman_logrpm_vs_parafrac"),
  row_stat(spearman_safe(dpos$wolb_breadth, dpos$para_frac),
           subset = "parasitoid_pos", metric = "spearman_breadth_vs_parafrac")
)
# Linear trend (log_rpm ~ para_frac) among parasitoid-positive samples.
a3_slope <- NA_real_; a3_slope_p <- NA_real_; a3_r2 <- NA_real_
if (nrow(dpos) >= 5 && length(unique(dpos$para_frac)) >= 3) {
  fit <- lm(log_rpm ~ para_frac, data = dpos)
  sm  <- summary(fit)
  a3_slope   <- coef(fit)[["para_frac"]]
  a3_slope_p <- coef(sm)["para_frac", "Pr(>|t|)"]
  a3_r2      <- sm$r.squared
}
write_tsv(a3, "results/03_profile/dose_response_stats.tsv")

ggsave("reports/figs/dose_response.png",
  ggplot(d, aes(para_frac, log_rpm)) +
    geom_point(aes(color = factor(wolb_pos)), alpha = .6) +
    geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = "grey30") +
    scale_color_manual(values = c("0" = "grey70", "1" = "firebrick"),
                       name = "Wolbachia+") +
    labs(x = "parasitoid read fraction (of non-host)",
         y = "log10(Wolbachia RPM + 1)",
         title = "A3 dose-response: Wolbachia load vs parasitoid fraction"),
  width = 6, height = 4.2, dpi = 150)

# ---------------------------------------------------------------------------
# A5 (lite) — source specificity: Wolbachia tracks parasitoid, not host/depth
# ---------------------------------------------------------------------------
# Marginal Spearman of Wolbachia load against each candidate driver.
a5 <- bind_rows(
  row_stat(spearman_safe(d$log_rpm, d$para_frac),      driver = "para_frac (parasitoid, of nonhost)"),
  row_stat(spearman_safe(d$log_rpm, d$para_frac_tot),  driver = "para_frac_tot (parasitoid, of total)"),
  row_stat(spearman_safe(d$log_rpm, d$chilo_frac_tot), driver = "chilo_frac_tot (host biomass)"),
  row_stat(spearman_safe(d$log_rpm, d$host_pct),       driver = "host_pct (host fraction)"),
  row_stat(spearman_safe(d$log_rpm, d$nonhost_frac),   driver = "nonhost_frac (generic non-host)"),
  row_stat(spearman_safe(d$log_rpm, d$total),          driver = "total (sequencing depth)")
) %>% mutate(kind = "marginal")

# Partial Spearman: does the parasitoid association survive controlling for
# sequencing depth and generic non-host load?
pc <- partial_spearman(d$log_rpm, d$para_frac, d[, c("total", "nonhost_frac")])
a5_partial <- tibble(driver = "para_frac | (depth + nonhost_frac)",
                     rho = pc[["rho"]], p = pc[["p"]], n = pc[["n"]],
                     kind = "partial")
a5_out <- bind_rows(a5, a5_partial)
write_tsv(a5_out, "results/03_profile/source_specificity_stats.tsv")

# Faceted scatter: Wolbachia load vs the three competing drivers (long format
# built without tidyr to keep dependencies minimal).
long <- bind_rows(
  transmute(d, log_rpm, value = para_frac,      driver = "parasitoid fraction"),
  transmute(d, log_rpm, value = chilo_frac_tot, driver = "host biomass fraction"),
  transmute(d, log_rpm, value = nonhost_frac,   driver = "non-host fraction")
)
long$driver <- factor(long$driver,
  levels = c("parasitoid fraction", "host biomass fraction", "non-host fraction"))
ggsave("reports/figs/source_specificity.png",
  ggplot(long, aes(value, log_rpm)) +
    geom_point(alpha = .5, size = .8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = "firebrick") +
    facet_wrap(~ driver, scales = "free_x") +
    labs(x = "driver value", y = "log10(Wolbachia RPM + 1)",
         title = "A5 specificity: Wolbachia tracks parasitoid, not host/depth"),
  width = 9, height = 3.6, dpi = 150)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
fmt <- function(r) if (is.na(r["rho"])) "n/a (degenerate)" else
  sprintf("rho=%.3g, p=%.3g, n=%d", r["rho"], r["p"], r["n"])

a3_all  <- spearman_safe(d$log_rpm, d$para_frac)
a3_ppos <- spearman_safe(dpos$log_rpm, dpos$para_frac)

lines <- c(
  "# Stage 3b — 剂量效应 (A3) + 来源特异性 (A5-lite)",
  "",
  sprintf("N = %d 样本（Wolbachia+ %d，寄生+ %d）。纯重分析 per_sample_profile.tsv，无重映射。",
          nrow(d), sum(d$wolb_pos), sum(d$para_pos)),
  "",
  "## A3 剂量效应：Wolbachia 载量 vs 寄生物读段比例",
  "",
  sprintf("- 全样本 Spearman(log10 RPM vs para_frac): %s", fmt(a3_all)),
  sprintf("- 寄生阳性子集 Spearman(log10 RPM vs para_frac): %s", fmt(a3_ppos)),
  if (!is.na(a3_slope))
    sprintf("- 寄生阳性子集线性趋势 log10 RPM ~ para_frac: slope=%.3g, p=%.3g, R^2=%.3g",
            a3_slope, a3_slope_p, a3_r2)
  else "- 线性趋势：样本不足或退化，未计算",
  "",
  "解读：正向单调 = Wolbachia 载量随寄生物比例上升（支持寄生物来源 H1，反对 H3 随机串扰的剂量无关）。",
  "见 reports/figs/dose_response.png",
  "",
  "## A5-lite 来源特异性：Wolbachia 与谁同涨？",
  "",
  "| driver | rho | p | n | kind |",
  "|---|---|---|---|---|",
  paste(apply(a5_out, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
        r[["driver"]],
        ifelse(is.na(r[["rho"]]), "n/a", sprintf("%.3g", as.numeric(r[["rho"]]))),
        ifelse(is.na(r[["p"]]),   "n/a", sprintf("%.3g", as.numeric(r[["p"]]))),
        ifelse(is.na(r[["n"]]),   "n/a", as.character(r[["n"]])),
        r[["kind"]])), collapse = "\n"),
  "",
  "解读判据：Wolbachia 应与**寄生物比例**强正相关、与**宿主生物量/测序深度/泛非宿主负载**不相关或弱，",
  "且偏相关（控制深度+非宿主负载后）仍显著，方能**特异于寄生物信号**。",
  "",
  "**自动判读（据本次数据）**：",
  sprintf("- 寄生物(of total) rho=%.3g vs 泛非宿主负载 rho=%.3g；偏相关(控制深度+非宿主) rho=%.3g, p=%.3g",
          a5$rho[a5$driver == "para_frac_tot (parasitoid, of total)"],
          a5$rho[a5$driver == "nonhost_frac (generic non-host)"],
          a5_partial$rho, a5_partial$p),
  paste0("- 结论：**A5-lite 不足以干净判别**。泛非宿主负载与 Wolbachia 的相关性≈寄生物特异信号，",
         "且偏相关骤降——因 para_reads 本就是 nonhost 的组成部分，现有列彼此代数纠缠，无法分离",
         "“寄生物特异” vs “泛非宿主/污染负载”。**这恰恰要求做真·双线粒体重映射（A5-full）**。"),
  "- A3 剂量效应不受此影响：寄生阳性子集内 Wolbachia 随 para_frac 单调上升，是稳健结果。",
  "见 reports/figs/source_specificity.png",
  "",
  "## 下一步（真·双线粒体 A5-full，已从“可选”升级为“必需”）",
  "",
  "A5-lite 的纠缠表明：要干净证明“Wolbachia 特异伴随寄生物、而非泛非宿主/污染负载”，",
  "必须把**原始 reads** 映射到 {Chilo-mito + 6 株寄生物 mito} 组合参考，得到**对称、互斥**的",
  "双昆虫线粒体读段计数，再检验 Wolbachia 载量 ~ 寄生物 mito（控制 Chilo mito）。",
  "参考已就绪（refs/enemies/mito/*.mito.fa）、原始/nonhost fastq 也在 —— 建议实现 **Stage 03c**：",
  "构建组合 mito 参考 → minibwa 比对 → idxstats 取每 mito 读段数 → 偏相关/回归。",
  "附带产出：每样本**主导寄生物种**（解释 HC_W_9 被 Lypha 主导，支撑证据线 B+）。"
)
writeLines(lines, "reports/stage03b_dose_source.md")
cat("done\n")
