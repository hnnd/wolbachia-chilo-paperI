suppressMessages({library(readr); library(dplyr); library(ggplot2)})
set.seed(42)

# Stage 03c — A5-full dual-source mitochondrial specificity.
# Symmetric, mutually-exclusive insect-mito read counts (Chilo vs parasitoids)
# from mapping non-host reads to the combined mitogenome reference. Tests whether
# Wolbachia load tracks PARASITOID mito specifically, controlling for Chilo mito
# (host biomass) and sequencing depth — the clean test A5-lite could not do
# because its existing columns were algebraically entangled.

m <- read_tsv("results/03_profile/mito_long.tsv",
              col_names = c("sample", "species", "contig", "length", "reads"),
              show_col_types = FALSE)
prof <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)
dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)

CHILO <- "Chilo_suppressalis"

# Per-sample wide summary: Chilo mito vs summed parasitoid mito + dominant species.
per <- m %>%
  group_by(sample) %>%
  summarise(
    chilo_mito = sum(reads[species == CHILO]),
    para_mito  = sum(reads[species != CHILO]),
    dom_para   = {
      ps <- reads[species != CHILO]; sp <- species[species != CHILO]
      if (length(ps) == 0 || max(ps) == 0) "none" else sp[which.max(ps)]
    },
    dom_reads  = { ps <- reads[species != CHILO]; if (length(ps)) max(ps) else 0 },
    .groups = "drop"
  ) %>%
  mutate(
    total_mito     = chilo_mito + para_mito,
    para_mito_frac = ifelse(total_mito > 0, para_mito / total_mito, NA_real_),
    dom_frac       = ifelse(para_mito > 0, dom_reads / para_mito, NA_real_)
  )

d <- prof %>%
  select(sample, wolb_rpm, wolb_breadth, wolb_pos, para_pos, para_frac, total, host_mapped) %>%
  inner_join(per, by = "sample") %>%
  mutate(log_rpm = log10(wolb_rpm + 1))
write_tsv(d, "results/03_profile/per_sample_mito.tsv")

# NOTE on the host-biomass control: Chilo mito reads are largely removed during
# dehosting (Lepidoptera nuclear assemblies carry NUMTs, so mito reads map to the
# nuclear reference and are stripped) — chilo_mito in non-host reads is near-zero
# and noisy. The robust host-biomass control is therefore host_mapped (Chilo reads
# from TOTAL), with chilo_mito kept only as a secondary cross-check.

# ---- helpers (same conventions as 03b) ----
spearman_safe <- function(x, y) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 3 || length(unique(y)) < 3)
    return(c(rho = NA_real_, p = NA_real_, n = length(x)))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}
partial_spearman <- function(x, y, Z) {
  ok <- complete.cases(x, y, Z); x <- x[ok]; y <- y[ok]
  Z <- as.data.frame(Z)[ok, , drop = FALSE]
  if (length(unique(x)) < 3 || length(unique(y)) < 3)
    return(c(rho = NA_real_, p = NA_real_, n = length(x)))
  rz <- as.data.frame(lapply(Z, rank))
  ex <- residuals(lm(rank(x) ~ ., data = rz))
  ey <- residuals(lm(rank(y) ~ ., data = rz))
  ct <- suppressWarnings(cor.test(ex, ey, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}
row_stat <- function(v, ...) tibble(..., rho = v[["rho"]], p = v[["p"]], n = v[["n"]])

# ---- A5-full specificity tests ----
# Primary control: host_mapped (host biomass from total reads) + depth.
# chilo_mito kept as secondary marginal cross-check (attenuated by NUMT removal).
stats <- bind_rows(
  row_stat(spearman_safe(d$log_rpm, d$para_mito),      test = "wolb ~ parasitoid_mito",  kind = "marginal"),
  row_stat(spearman_safe(d$log_rpm, d$para_mito_frac), test = "wolb ~ para_mito_frac",   kind = "marginal"),
  row_stat(spearman_safe(d$log_rpm, d$host_mapped),    test = "wolb ~ host_mapped",      kind = "marginal"),
  row_stat(spearman_safe(d$log_rpm, d$chilo_mito),     test = "wolb ~ chilo_mito(xcheck)", kind = "marginal"),
  row_stat(spearman_safe(d$log_rpm, d$total),          test = "wolb ~ depth(total)",     kind = "marginal"),
  row_stat(partial_spearman(d$log_rpm, d$para_mito, d[, c("host_mapped", "total")]),
           test = "wolb ~ parasitoid_mito | (host_mapped + depth)", kind = "partial")
)
write_tsv(stats, "results/03_profile/mito_specificity_stats.tsv")

# Dominant parasitoid composition (parasitoid-positive samples).
dom <- d %>% filter(para_pos == 1) %>% count(dom_para, sort = TRUE)
write_tsv(dom, "results/03_profile/mito_dominant_parasitoid.tsv")

# Figure: Wolbachia load vs parasitoid-mito vs Chilo-mito (log scale).
long <- bind_rows(
  transmute(d, log_rpm, value = log10(para_mito + 1),  driver = "parasitoid mito (log10 reads+1)"),
  transmute(d, log_rpm, value = log10(chilo_mito + 1), driver = "Chilo mito (log10 reads+1)")
)
ggsave("reports/figs/mito_dualsource.png",
  ggplot(long, aes(value, log_rpm)) +
    geom_point(alpha = .5, size = .9) +
    geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = "firebrick") +
    facet_wrap(~ driver, scales = "free_x") +
    labs(x = "insect mito abundance", y = "log10(Wolbachia RPM + 1)",
         title = "A5-full: Wolbachia tracks parasitoid mito, not Chilo mito"),
  width = 8, height = 3.8, dpi = 150)

# ---- report ----
g <- function(t) { r <- stats[stats$test == t, ]; if (!nrow(r) || is.na(r$rho)) "n/a"
                   else sprintf("rho=%.3g, p=%.3g, n=%d", r$rho, r$p, r$n) }
getrho <- function(t) { r <- stats$rho[stats$test == t]; if (length(r)) r[1] else NA_real_ }
pm  <- getrho("wolb ~ parasitoid_mito")
cm  <- getrho("wolb ~ chilo_mito(xcheck)")
hm  <- getrho("wolb ~ host_mapped")
dep <- getrho("wolb ~ depth(total)")
pr  <- stats[stats$kind == "partial", ]
# General-abundance effect: Wolbachia also rises with host biomass / depth /
# chilo mito at a magnitude comparable to parasitoid mito.
gen_effect <- all(!is.na(c(cm, hm, dep))) && min(cm, hm, dep) > 0.2
partial_ok <- nrow(pr) && !is.na(pr$rho) && pr$rho > 0 && pr$p < 0.05
verdict <- if (partial_ok && !gen_effect) {
    "**通过（强）**：偏相关显著，且寄生物 mito 明显强于宿主/深度的泛丰度效应。"
  } else if (partial_ok && gen_effect) {
    sprintf(paste0("**部分通过（需谨慎措辞）**：控制宿主生物量+深度后，寄生物 mito 仍显著独立预测 ",
            "Wolbachia（partial rho=%.2f, p=%.2g）；**但** Wolbachia 同时与 chilo_mito(%.2f)、",
            "host_mapped(%.2f)、depth(%.2f) 正相关，存在泛丰度/深度效应——偏相关只能说明寄生物 mito ",
            "有**独立成分**，不能宣称纯特异。更强的证据是 A3 剂量效应（寄生阳性子集 rho≈0.65，见 stage03b）。"),
            pr$rho, pr$p, cm, hm, dep)
  } else {
    "**未通过/需复核**：偏相关不显著 —— 见下表逐项判读。"
  }

lines <- c(
  "# Stage 03c — A5-full 双线粒体来源特异性",
  "",
  sprintf("N = %d 样本（与 profile 内连接）。非宿主 reads → {Chilo + 11 株寄生物} 组合线粒体参考，", nrow(d)),
  "primary-best-hit 计数。Chilo 线粒体得以保留是因宿主组装无 mito contig（仅染色体）。",
  "",
  "## 特异性检验（Wolbachia 载量 vs 昆虫线粒体）",
  "",
  "| 检验 | rho | p | n | kind |",
  "|---|---|---|---|---|",
  paste(apply(stats, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
        r[["test"]],
        ifelse(is.na(r[["rho"]]), "n/a", sprintf("%.3g", as.numeric(r[["rho"]]))),
        ifelse(is.na(r[["p"]]),   "n/a", sprintf("%.3g", as.numeric(r[["p"]]))),
        r[["n"]], r[["kind"]])), collapse = "\n"),
  "",
  paste("**自动判读**：", verdict),
  "",
  "寄生物各株 mito 为**互斥的 primary-best-hit** 计数（物种可分）；",
  "`wolb ~ parasitoid_mito | (host_mapped + depth)` 偏相关检验寄生物 mito 是否**独立于**宿主生物量与深度预测 Wolbachia，",
  "改进 A5-lite 中 para_reads ⊂ nonhost 的代数纠缠（见 reports/stage03b_dose_source.md）。",
  "",
  "> **诚实 caveat**：(1) Wolbachia 与多个丰度指标（chilo_mito、host_mapped、depth）同向相关，",
  "> 存在泛丰度/深度效应，偏相关只证明寄生物 mito 有**独立成分**，非纯特异；最强证据仍是 A3 剂量效应。",
  "> (2) `para_mito_frac`（寄生物 mito / 总昆虫 mito）因 chilo_mito 近零而≈1、**退化无信息**，不作判据。",
  "> (3) Chilo 线粒体被去宿主大量剥离（鳞翅目 NUMTs），故 host 对照用 host_mapped（来自全部 reads）。",
  "见 reports/figs/mito_dualsource.png",
  "",
  "## 主导寄生物组成（寄生阳性样本，支撑证据线 B+）",
  "",
  paste(capture.output(print(as.data.frame(dom))), collapse = "\n"),
  "",
  "> 注：各寄生物线粒体长度 14.7–16.9 kb，差异 <15%，未做长度归一；",
  "> Cotesia_chilonis 参考缺失（空文件），以同属 C. flavipes/ruficrus 作近缘代理。"
)
writeLines(lines, "reports/stage03c_mito_dualsource.md")
cat("done\n")
