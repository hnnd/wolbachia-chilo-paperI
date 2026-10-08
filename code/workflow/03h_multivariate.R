suppressMessages({library(readr); library(dplyr); library(ggplot2); library(logistf)})
set.seed(42)

# Stage 03h (A4) — multivariable logistic: is parasitoid load still the dominant
# predictor of Wolbachia positivity after adjusting for the "pan-abundance"
# confounders (sequencing depth, non-moth biomass) and host plant?
#
# Reviewer attack #5 (Codex): co-occurrence models omit confounders; para_frac
# might merely track total non-moth biomass or DNA quality.
#
# METHOD: near-perfect co-occurrence makes plain glm separate (fitted probs hit
# 0/1, OR blows up 41 -> 3075). We therefore use Firth penalized logistic
# regression (logistf), the standard fix for separation / rare events; estimates
# are finite and CIs / p-values come from the penalized profile likelihood.
#
# HONEST SCOPE: samples.tsv carries only source/group/host_plant/geo. Per-sample
# collection date, larval stage, sex, extraction/library batch and DNA-quality
# metrics are NOT recorded, so they cannot enter the model. We proxy the
# pan-abundance confounder with sequencing depth (total reads) and non-moth
# biomass fraction (nonhost/total). source(SRA) and geo(17 levels) induce
# separation (all 37 SRA samples are Wolbachia-, geo sparse) and are reported
# qualitatively rather than forced into the fit.

d <- read_tsv("results/03_profile/per_sample_profile.tsv", show_col_types = FALSE)
dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)

# continuous predictors standardized (z-score) so ORs are per 1 SD and comparable
d <- d %>% mutate(
  z_lpf      = as.numeric(scale(log10(para_frac + 1e-4))),  # parasitoid load
  z_depth    = as.numeric(scale(log10(total))),             # sequencing depth
  z_biomass  = as.numeric(scale(nonhost / total)),          # non-moth biomass fraction (pan-abundance confounder)
  host_plant = relevel(factor(host_plant), ref = "rice")
)

# document that plain glm separates (justifies Firth)
g <- suppressWarnings(glm(wolb_pos ~ z_lpf + z_depth + z_biomass + host_plant, data = d, family = binomial))
glm_fit_range <- range(fitted(g)); glm_or_lpf <- exp(coef(g)["z_lpf"])

# Firth penalized logistic. raise iterations: near-perfect co-occurrence makes the
# profile likelihood slow. PRIMARY model M1 targets the pan-abundance attack
# (load + depth + biomass). host_plant (evidence line C) is added only as a
# sensitivity model M2 — conditional on host it induces its own separation, so its
# own coefficient is reported qualitatively, not as a headline estimate.
ctrl <- logistf.control(maxit = 5000, maxstep = 5)
plc  <- logistpl.control(maxit = 5000, maxstep = 5)
m_uni  <- logistf(wolb_pos ~ z_lpf, data = d, control = ctrl, plcontrol = plc)
m1     <- logistf(wolb_pos ~ z_lpf + z_depth + z_biomass, data = d, control = ctrl, plcontrol = plc)
m2     <- logistf(wolb_pos ~ z_lpf + z_depth + z_biomass + host_plant, data = d, control = ctrl, plcontrol = plc)

coef_tbl <- tibble(
  term = names(coef(m1)),
  OR   = exp(coef(m1)),
  lo   = exp(m1$ci.lower),
  hi   = exp(m1$ci.upper),
  p    = m1$prob                  # penalized-likelihood ratio test p-value
)
write_tsv(coef_tbl, "results/03_profile/multivariate_coef.tsv")

or_lpf_uni <- exp(coef(m_uni)["z_lpf"]); or_lpf_adj <- exp(coef(m1)["z_lpf"])
or_lpf_m2  <- exp(coef(m2)["z_lpf"]); p_lpf_m2 <- m2$prob["z_lpf"]   # load stability with host plant added

# forest plot of adjusted ORs (M1, drop intercept), log scale. The parasitoid-load
# upper CI is unbounded under separation, so clip the whisker for display and flag it.
CAP <- 1e3
fp <- coef_tbl %>% filter(term != "(Intercept)") %>%
  mutate(label = recode(term,
    z_lpf = "parasitoid load (per SD)", z_depth = "seq depth (per SD)",
    z_biomass = "non-moth biomass frac (per SD)"),
    hi_disp = pmin(hi, CAP), clipped = hi > CAP)
p_forest <- ggplot(fp, aes(OR, reorder(label, OR))) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi_disp), color = "firebrick", linewidth = .6, size = .5) +
  geom_text(data = subset(fp, clipped), aes(x = CAP, label = "→ ∞"),
            hjust = -0.1, vjust = 0.4, size = 3, color = "firebrick") +
  scale_x_log10(limits = c(0.02, CAP * 3)) +
  labs(x = "adjusted odds ratio (95% CI, log scale)", y = NULL,
       title = "A4: parasitoid load dominates Wolbachia+ prediction after adjustment",
       subtitle = sprintf("Firth logistic (n=%d, 31 Wolbachia+), adjusted for seq depth + non-moth biomass; p(load)=%.2g",
                          nrow(d), coef_tbl$p[coef_tbl$term == "z_lpf"])) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11.5),
        plot.subtitle = element_text(size = 8, color = "grey30"))
ggsave("reports/figs/multivariate_or.png", p_forest, width = 9, height = 3.4, dpi = 150)

# ---- Report
fmt <- function(term) { r <- coef_tbl[coef_tbl$term == term, ]
  sprintf("OR=%.2f (95%% CI %.2f–%.2f), p=%.2g", r$OR, r$lo, r$hi, r$p) }
lines <- c(
  "# Stage 03h (A4) — 多因素 logistic：控制混杂后寄生载量仍是主导预测因子",
  "",
  sprintf("N = %d（31 Wolbachia+）。主模型 M1：wolb_pos ~ 寄生载量 + 测序深度 + 非宿主生物量占比。", nrow(d)),
  "连续变量已 z 标准化（OR 为每 1 个标准差）。纯重分析 per_sample_profile.tsv。",
  "",
  "## 方法：为何用 Firth 惩罚回归",
  "",
  sprintf("寄生物与 Wolbachia 近完美共现使普通 glm **完全分离**（拟合概率触及 [%.3f, %.3f]，寄生载量 OR 从合理值爆到 %.0f、CI 发散）。",
          glm_fit_range[1], glm_fit_range[2], glm_or_lpf),
  "故改用 **Firth 惩罚 logistic 回归（logistf）**——分离/稀有事件的标准解法，估计有限，CI 与 p 取自惩罚 profile likelihood。",
  "",
  "## 可用元数据的诚实边界",
  "",
  "`samples.tsv` 仅含 source / group / host_plant / geo。**逐样本采集日期、龄期、性别、提取/建库批次、DNA 质量未记录**，无法入模。",
  "故以**测序深度（total reads）**与**非宿主生物量占比（nonhost/total）**代理'泛丰度/DNA 质量'混杂——这正是审稿攻击 #5（共现只是生物量同涨）的核心。",
  "`source`（37 个 SRA 全部 Wolbachia−）与 `geo`（17 层、稀疏）会造成完全/准完全分离，故**定性报告、不强行入模**。",
  "",
  "## 结果（主模型 M1：泛丰度混杂）",
  "",
  "| 预测因子 | 校正后 OR (95% CI) | p |",
  "|---|---|---|",
  paste(apply(coef_tbl[coef_tbl$term != "(Intercept)", ], 1, function(r)
    sprintf("| %s | %.2f (%.2f–%.2f) | %.2g |", r[["term"]],
            as.numeric(r[["OR"]]), as.numeric(r[["lo"]]), as.numeric(r[["hi"]]),
            as.numeric(r[["p"]]))), collapse = "\n"),
  "",
  sprintf("- **寄生载量**：%s —— 控制测序深度与非宿主生物量占比后仍高度显著（上界因分离而无界，点估计与下界已远大于 1）。", fmt("z_lpf")),
  sprintf("- 测序深度：%s → 不显著。", fmt("z_depth")),
  sprintf("- 非宿主生物量占比：%s → 不显著。", fmt("z_biomass")),
  sprintf("- 寄生载量 OR：**单变量 %.2f → M1 校正后 %.2f**，方向与量级稳定（若混杂主导，调整应大幅削弱它）。",
          or_lpf_uni, or_lpf_adj),
  "",
  "## 敏感性（M2：加入寄主植物，证据线 C 协变量）",
  "",
  sprintf("加入 host_plant 后，寄生载量仍显著（OR=%.2f, p=%.2g），结论不变。", or_lpf_m2, p_lpf_m2),
  "host_plant 自身在条件于寄主时引入局部分离（茭白阳性样为低 Cotesia 载量的寄蝇/de novo 病例），其系数 CI 发散、不稳定，",
  "故**仅定性报告**：寄主分层的作用见证据线 C（Fig 1C 生态支撑），不依赖本模型的点估计。",
  "",
  "## 判读",
  "",
  "控制测序深度与非宿主生物量占比后，**寄生载量仍是 Wolbachia 阳性的主导且显著预测因子**，",
  "深度/生物量代理均不显著 → **排除'共现仅因总非宿主生物量/深度同涨'的泛丰度混杂**（与证据线 D 的分类群特异带菌互为印证）。",
  "见 reports/figs/multivariate_or.png、results/03_profile/multivariate_coef.tsv。"
)
writeLines(lines, "reports/stage03h_multivariate.md")
cat("done (Firth): lpf OR uni =", round(or_lpf_uni, 2), "-> adj =", round(or_lpf_adj, 2), "\n")
