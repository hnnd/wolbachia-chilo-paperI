#!/usr/bin/env Rscript
# Fig 4 — Evidence line D: parasitoid taxon-specific Wolbachia carriage.
# The anti-confounder result: Wolbachia carriage tracks parasitoid LINEAGE, not
# biomass. Ichneumonidae carry the largest parasitoid biomass yet zero Wolbachia;
# Braconidae (Cotesia) and Tachinidae carry it. Source: Stage 03f de novo + remote nt.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })

figdir <- "reports/figs"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
d <- read.delim("results/06_denovo/fig4_carriage.tsv", header = FALSE,
                col.names = c("sample","family","total_para","wolb_rpm","wolb_pos"))
d$total_para <- as.numeric(d$total_para); d$wolb_rpm <- as.numeric(d$wolb_rpm)
fam_lv <- c("Braconidae (Cotesia)","Tachinidae (fly)","Ichneumonidae")
d$family <- factor(d$family, levels = fam_lv)
pal <- setNames(c("#0072B2","#E69F00","#999999"), fam_lv)

# ---- A: biomass vs Wolbachia load, colored by family ----------------------
d$rpm_floor <- pmax(d$wolb_rpm, 0.1)   # pseudocount for log axis
pA <- ggplot(d, aes(total_para, rpm_floor, fill = family)) +
  annotate("rect", xmin = 5e4, xmax = 3e5, ymin = 0.07, ymax = 5,
           fill = "grey80", alpha = .35) +
  annotate("text", x = 1.2e5, y = 0.13, label = "high biomass,\nzero Wolbachia",
           size = 3, fontface = "italic", color = "grey25") +
  geom_hline(yintercept = 200, linetype = 2, color = "grey60") +
  annotate("text", x = 1.5e3, y = 260, label = "positivity threshold (rpm 200)",
           size = 2.8, color = "grey45", hjust = 0) +
  geom_point(size = 4, shape = 21, color = "grey20", alpha = .9) +
  scale_x_log10(labels = scales::comma) + annotation_logticks(sides = "b") +
  scale_y_log10() +
  scale_fill_manual(values = pal, name = "Parasitoid family") +
  labs(x = "parasitoid biomass (assembled-mito reads)",
       y = "Wolbachia load (rpm, log)",
       title = "Carriage tracks parasitoid lineage, not biomass",
       subtitle = "Ichneumonidae: 810k reads total yet 0/11 Wolbachia+; Cotesia/Tachinidae carry it") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"),
        legend.position = c(0.02,0.98), legend.justification = c(0,1),
        legend.background = element_rect(fill = alpha("white",.75), color=NA),
        legend.key.size = unit(0.4,"cm"), legend.title = element_text(size=8.5),
        legend.text = element_text(size=8))

# ---- B: carriage rate by family (with n and total biomass) ----------------
s <- d %>% group_by(family) %>%
  summarise(n = n(), pos = sum(wolb_pos), rate = pos/n,
            biomass = sum(total_para), .groups = "drop")
pB <- ggplot(s, aes(family, rate, fill = family)) +
  geom_col(width = .65, color = "grey20") +
  geom_text(aes(label = sprintf("%d/%d", pos, n)), vjust = -0.4, size = 3.4) +
  geom_text(aes(y = -0.06, label = sprintf("biomass\n%s reads", scales::comma(biomass))),
            size = 2.7, color = "grey30", lineheight = .9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(-0.12, 1.1), breaks = c(0,.5,1),
                     labels = scales::percent) +
  labs(x = NULL, y = "Wolbachia+ rate",
       title = "Wolbachia carriage by parasitoid family",
       subtitle = "highest-biomass lineage (Ichneumonidae) is Wolbachia-free") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"),
        axis.text.x = element_text(size = 8.5, face = c("bold","bold","bold")))

fig <- (pA | pB) +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure 4. Parasitoid taxon-specific Wolbachia carriage rules out a biomass/depth confounder",
    theme = theme(plot.title = element_text(face="bold", size=12))) &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(figdir, "Fig3_taxon_carriage.png"), fig, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(figdir, "Fig3_taxon_carriage.pdf"), fig, width = 12, height = 5.5)
cat("[Fig3] wrote", file.path(figdir, "Fig3_taxon_carriage.png/.pdf"), "\n")
