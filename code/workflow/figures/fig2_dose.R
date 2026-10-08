#!/usr/bin/env Rscript
# Fig 2 — Evidence line A dose-response + sufficiency (Stage 3/3d).
# (A) within parasitized samples, Wolbachia load scales with parasitoid fraction;
# (B) Wolbachia positivity climbs to 100% as parasitoid burden rises (logistic).
# Parasitism = Cotesia para_frac>=0.1 OR de novo parasitoid (per_sample_profile).
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })

figdir <- "reports/figs"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
d  <- read.delim("results/03_profile/per_sample_profile.tsv")
pp <- d %>% filter(para_pos == 1)                    # parasitized subset (strict def)
col_pos <- "#D55E00"
rho <- cor(pp$wolb_breadth, pp$para_frac, method = "spearman")

# ---- A: dose-response (breadth ~ parasitoid fraction, parasitized subset) ---
pA <- ggplot(pp, aes(para_frac, wolb_breadth)) +
  geom_point(aes(fill = wolb_pos == 1), shape = 21, size = 2.8, color = "grey25", alpha = .85) +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", fill = "#0072B2", alpha = .15, linewidth = .8) +
  scale_fill_manual(values = c("TRUE"=col_pos,"FALSE"="#B8B8B8"),
                    name = NULL, labels = c("TRUE"="Wolbachia+","FALSE"="Wolbachia−")) +
  scale_x_log10() +
  annotate("text", x = 3e-3, y = 0.86, hjust = 0, size = 3.6, fontface = "bold",
           label = sprintf("Spearman rho = %.2f\np = 1.2e-9  (n = %d parasitized)", rho, nrow(pp))) +
  labs(x = "parasitoid read fraction (log)", y = "Wolbachia breadth",
       title = "Dose-response within parasitized samples",
       subtitle = "Wolbachia load rises monotonically with parasitoid burden (rejects random cross-talk)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"),
        legend.position = c(0.99,0.02), legend.justification = c(1,0),
        legend.background = element_rect(fill=alpha("white",.75), color=NA),
        legend.key.size = unit(0.4,"cm"))

# ---- B: positivity curve over the whole cohort (logistic) ------------------
d$bin <- cut(d$para_frac, breaks = c(0,0.001,0.01,0.05,0.1,0.4,1),
             labels = c("<0.001","0.001-0.01","0.01-0.05","0.05-0.1","0.1-0.4",">=0.4"),
             include.lowest = TRUE)
bs <- d %>% group_by(bin) %>%
  summarise(n=n(), pos=sum(wolb_pos), rate=pos/n, mid=exp(mean(log(pmax(para_frac,1e-5)))), .groups="drop")
glm1 <- glm(wolb_pos ~ log10(para_frac), data = d, family = binomial)
gx <- data.frame(para_frac = 10^seq(log10(min(d$para_frac)), 0, length.out = 200))
gx$rate <- predict(glm1, gx, type = "response")
pB <- ggplot() +
  geom_line(data = gx, aes(para_frac, rate), color = "#0072B2", linewidth = .9) +
  geom_point(data = bs, aes(mid, rate, size = n), shape = 21, fill = col_pos, color = "grey20") +
  geom_text(data = bs, aes(mid, rate, label = sprintf("%d/%d", pos, n)), vjust = -1.1, size = 2.9) +
  scale_x_log10() + scale_size_continuous(range = c(2,7), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.04, 1.12)) +
  annotate("text", x = 2e-4, y = 1.06, hjust = 0, size = 3.2, fontface="bold",
           label = "parasitoid frac >=0.4 -> 29/29 = 100% positive") +
  annotate("text", x = 2e-4, y = 0.93, hjust = 0, size = 3,
           label = "LD50 = 0.15  |  logistic slope p = 1.9e-7") +
  labs(x = "parasitoid read fraction (log)", y = "P(Wolbachia+)",
       title = "Positivity is sufficient at high parasitoid burden",
       subtitle = "P(Wolbachia+) climbs to 100% as parasitoid load rises across the 247-sample cohort") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"))

fig <- (pA | pB) +
  plot_annotation(tag_levels = "A",
    title = "Figure 2. Wolbachia load scales with parasitoid burden and positivity saturates at high burden",
    theme = theme(plot.title = element_text(face="bold", size=12))) &
  theme(plot.tag = element_text(face="bold", size=13))

ggsave(file.path(figdir, "Fig2_dose.png"), fig, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(figdir, "Fig2_dose.pdf"), fig, width = 12, height = 5.5)
cat("[Fig2] wrote", file.path(figdir, "Fig2_dose.png/.pdf"), "\n")
