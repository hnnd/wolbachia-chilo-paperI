#!/usr/bin/env Rscript
# Fig 1 — Evidence line A: Wolbachia co-occurs only with parasitoid signal.
# (A) per-sample scatter with the empty "Wolbachia+ / parasitoid-" quadrant,
# (B) 2x2 contingency (Fisher), (C) host-plant stratification. Source: Stage 3.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })

figdir <- "reports/figs"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
d <- read.delim("results/03_profile/per_sample_profile.tsv")
# de novo parasitoid flag (Stage 03f): rescues fly/ichneumonid-parasitized samples
# whose Cotesia-referenced para_frac is low. para_pos in the profile already folds
# this in; here we recover the flag to MARK those samples in the scatter.
dn <- tryCatch(read.delim("results/06_denovo/per_sample_denovo_id.tsv"), error=function(e) NULL)
para_genera <- c("Cotesia","Pexopsis","Scambus","Neotrichoporoides","Lypha","Brachymeria",
                 "Trichogramma","Trichomalopsis","Pachyneuron","Tetrastichus")
d$denovo_para <- 0
if (!is.null(dn)) { dp <- dn$sample[dn$denovo_genus %in% para_genera & dn$mito_len >= 1000]
                    d$denovo_para <- as.integer(d$sample %in% dp) }
d$wolb <- factor(ifelse(d$wolb_pos == 1, "Wolbachia+", "Wolbachia−"),
                 levels = c("Wolbachia+","Wolbachia−"))
# parasitism route for point shape: Cotesia-threshold vs de-novo-rescued vs none
d$route <- factor(ifelse(d$para_frac >= 0.1, "Cotesia >=10%",
                  ifelse(d$denovo_para == 1, "de novo parasitoid", "not parasitized")),
                  levels = c("Cotesia >=10%","de novo parasitoid","not parasitized"))
pal <- c("Wolbachia+" = "#D55E00", "Wolbachia−" = "#B8B8B8")
BREADTH_THR <- 0.05; PARA_THR <- 0.1   # breadth is the discriminating axis (rpm alone false-positives, e.g. SAMN10449037)

# ---- A: para fraction vs Wolbachia genome breadth -------------------------
# breadth cleanly separates: positives >=0.089, negatives <=0.030 (0 overlap),
# so every Wolbachia+ point sits above the 5% line and every Wolbachia- below it.
hc <- d[d$sample == "HC_W_9", ]
pA <- ggplot(d, aes(para_frac + 1e-5, wolb_breadth)) +
  annotate("rect", xmin = 1e-5, xmax = PARA_THR, ymin = 0, ymax = Inf,
           fill = "grey85", alpha = .4) +
  geom_hline(yintercept = BREADTH_THR, linetype = 2, color = "grey55") +
  annotate("text", x = 5e-4, y = BREADTH_THR, vjust = -0.5, hjust = 0, size = 2.6,
           color = "grey45", label = "positivity threshold (5% breadth)", fontface = "italic") +
  geom_vline(xintercept = PARA_THR, linetype = 3, color = "grey55") +
  geom_point(aes(fill = wolb, shape = route), size = 2.7, color = "grey25", alpha = .85) +
  geom_point(data = hc, shape = 21, size = 5, stroke = 1.1, color = "#C0392B", fill = NA) +
  annotate("text", x = hc$para_frac, y = hc$wolb_breadth + 0.09, label = "HC_W_9\nfly-parasitized (de novo)",
           size = 2.7, color = "#C0392B", lineheight = .9) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.95)) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c("Cotesia >=10%"=21,"de novo parasitoid"=24,"not parasitized"=22),
                     name = "parasitism") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "parasitoid read fraction (Cotesia ref, log)", y = "Wolbachia genome breadth (% covered)",
       title = "Wolbachia appears only with parasitism",
       subtitle = "31/31 Wolbachia+ are parasitized (Cotesia >=10% or de novo parasitoid); empty Wolbachia+/unparasitized quadrant") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=7.6, color="grey30"),
        legend.position = c(0.01,0.99), legend.justification = c(0,1),
        legend.background = element_rect(fill=alpha("white",.75), color=NA),
        legend.key.size = unit(0.35,"cm"), legend.text = element_text(size=7.5),
        legend.title = element_text(size=8), legend.spacing.y = unit(0.02,"cm"))

# ---- B: 2x2 contingency ----------------------------------------------------
ct <- as.data.frame(table(Wolbachia = d$wolb,
                          Parasitoid = factor(ifelse(d$para_pos==1,"parasitoid+","parasitoid−"),
                                              levels=c("parasitoid+","parasitoid−"))))
ct$hl <- with(ct, Wolbachia=="Wolbachia+" & Parasitoid=="parasitoid−")
pB <- ggplot(ct, aes(Parasitoid, Wolbachia, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = Freq, color = hl), size = 7, fontface = "bold") +
  scale_fill_gradient(low = "#EAF2F8", high = "#2E86C1", guide = "none") +
  scale_color_manual(values = c("FALSE"="grey15","TRUE"="#C0392B"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Contingency: Wolbachia x parasitism",
       subtitle = "Fisher OR = Inf, p = 6.6e-29 (0 counter-examples; parasitism strictly defined)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 9.5, face = "bold"))

# ---- C: host-plant stratification -----------------------------------------
hs <- d %>% group_by(host_plant) %>%
  summarise(n = n(), pos = sum(wolb_pos), rate = pos/n, .groups="drop") %>%
  mutate(host = recode(host_plant, rice = "rice (_R)", wateroat = "water-oat (_W)"))
pC <- ggplot(hs, aes(host, rate)) +
  geom_col(width = .6, fill = "#D55E00", color = "grey20", alpha = .9) +
  geom_text(aes(label = sprintf("%d/%d\n%.1f%%", pos, n, 100*rate)), vjust = -0.25, size = 3.3, lineheight=.9) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.22)) +
  labs(x = NULL, y = "Wolbachia+ rate",
       title = "Host-plant stratification",
       subtitle = "rice-borne higher (ecological support, not source proof)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face="bold", size=11),
        plot.subtitle = element_text(size=8, color="grey30"),
        axis.text.x = element_text(size = 9.5))

fig <- pA + (pB / pC) + plot_layout(widths = c(1.3, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure 1. Wolbachia signal in field Chilo suppressalis co-occurs strictly with parasitoid signal",
    theme = theme(plot.title = element_text(face="bold", size=12))) &
  theme(plot.tag = element_text(face="bold", size=13))

ggsave(file.path(figdir, "Fig1_cooccurrence.png"), fig, width = 12, height = 6.5, dpi = 300)
ggsave(file.path(figdir, "Fig1_cooccurrence.pdf"), fig, width = 12, height = 6.5)
cat("[Fig1] wrote", file.path(figdir, "Fig1_cooccurrence.png/.pdf"), "\n")
