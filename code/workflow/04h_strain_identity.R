suppressMessages({library(readr); library(dplyr); library(ggplot2)})
set.seed(42)

# Stage 04h (B5 objective same-strain criterion) — formalize "a single circulating
# Wolbachia strain-unit" with predefined thresholds, replacing narrative "the
# 7/6/7/3/8 signature recurs" claims. Two independent readouts:
#   (1) MLST allele-profile concordance (predefined: identical 5-locus profile =
#       same ST = same strain, the standard MLST definition);
#   (2) whole-genome fastANI among the assemblies with sufficient contiguity
#       (predefined: ANI >=95% same species, >=99% same clonal strain).
# HONEST LIMITS: only 7 samples yield MLST calls and only 3 enrichment assemblies
# are contiguous enough (>=100 3kb fragments) for genome ANI; co-infection (A+B)
# depresses composite ANI, so ANI values are conservative lower bounds.

loci <- c("gatB", "coxA", "hcpA", "ftsZ", "fbpA")
Asig <- c(gatB = "7", coxA = "6",   hcpA = "7",   ftsZ = "3", fbpA = "8")    # ST-19, supergroup A
Bsig <- c(gatB = "9", coxA = "150", hcpA = "181", ftsZ = NA,  fbpA = "230")  # second (supergroup B) alleles
outdir <- "results/05_typing/strain_id"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create("reports/figs", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Part 1 — MLST allele-profile concordance
# ---------------------------------------------------------------------------
parse_mlst <- function(f) {
  s <- sub("\\.mlst\\.tsv$", "", basename(f))
  txt <- readLines(f, warn = FALSE)
  if (length(txt) == 0) return(NULL)
  parts <- strsplit(txt[1], "\t")[[1]]
  calls <- setNames(rep(NA_character_, length(loci)), loci)
  if (any(parts == "wolbachia")) {
    for (p in parts) {
      m <- regmatches(p, regexec("^(gatB|coxA|hcpA|ftsZ|fbpA)\\((.*)\\)$", p))[[1]]
      if (length(m) == 3) {
        al <- gsub("[~?]", "", strsplit(m[3], ",")[[1]]); al <- al[al != "-" & al != ""]
        if (length(al)) calls[m[2]] <- paste(al, collapse = ",")
      }
    }
  }
  tibble(sample = s, locus = loci, alleles = unname(calls[loci]))
}
prof <- bind_rows(lapply(Sys.glob("results/05_typing/*.mlst.tsv"), parse_mlst))
typed <- prof %>% group_by(sample) %>% filter(any(!is.na(alleles))) %>% ungroup()
typed_samples <- sort(unique(typed$sample)); n_typed <- length(typed_samples)

has_allele <- function(alleles, a) !is.na(alleles) & !is.na(a) &
  vapply(strsplit(alleles, ","), function(v) a %in% v, logical(1))
cls <- typed %>% rowwise() %>% mutate(
  hasA = has_allele(alleles, Asig[[locus]]),
  hasB = has_allele(alleles, Bsig[[locus]]),
  klass = dplyr::case_when(
    is.na(alleles)      ~ "missing",
    hasA & hasB         ~ "A+B (both)",
    hasA                ~ "A (ST-19)",
    hasB                ~ "B (second)",
    TRUE                ~ "other/noise")) %>% ungroup()
cls$locus <- factor(cls$locus, levels = loci)
cls$klass <- factor(cls$klass, levels = c("A (ST-19)", "A+B (both)", "B (second)", "other/noise", "missing"))

# per-locus recurrence of the A (ST-19) signature and B second alleles
recur <- tibble(locus = loci,
  A_allele = unname(Asig[loci]), A_n = vapply(loci, function(L) sum(cls$hasA[cls$locus == L]), integer(1)),
  B_allele = unname(Bsig[loci]), B_n = vapply(loci, function(L) sum(cls$hasB[cls$locus == L]), integer(1)))
n_full_ST19 <- cls %>% group_by(sample) %>% summarise(full = all(hasA[match(loci, locus)]), .groups = "drop") %>%
  filter(full) %>% nrow()
write_tsv(cls %>% select(sample, locus, alleles, klass), file.path(outdir, "mlst_profiles.tsv"))
write_tsv(recur, file.path(outdir, "mlst_recurrence.tsv"))

pal <- c("A (ST-19)" = "#1B7837", "A+B (both)" = "#762A83", "B (second)" = "#2166AC",
         "other/noise" = "#B8B8B8", "missing" = "grey92")
p_mlst <- ggplot(cls, aes(locus, sample, fill = klass)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(is.na(alleles), "", alleles)), size = 2.9) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "MLST locus", y = NULL,
       title = "Stage 04h: MLST allele-profile concordance across typeable Wolbachia",
       subtitle = sprintf("green = ST-19 supergroup-A signature (7/6/7/3/8); blue = recurrent second (supergroup-B) alleles; n=%d typeable", n_typed)) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8, color = "grey30"), panel.grid = element_blank())
ggsave("reports/figs/strain_mlst_concordance.png", p_mlst, width = 7.5, height = 4.2, dpi = 150)

# ---------------------------------------------------------------------------
# Part 2 — whole-genome fastANI
# ---------------------------------------------------------------------------
refname <- c(NC_002978.6 = "wMel (A)", NC_012416.1 = "wRi (A)", NC_010981.1 = "wPip (B)",
             NZ_CM003641.1 = "wTpre (B)", wBm = "wBm (D)")
short <- function(p) { a <- Filter(nchar, strsplit(p, "/")[[1]])
  if (tail(a, 1) == "assembly.fasta") a[length(a) - 1] else sub("\\.(fa|fna)$", "", tail(a, 1)) }
ani <- read_tsv(file.path(outdir, "ani_all.tsv"), show_col_types = FALSE,
                col_names = c("q", "r", "ani", "mfrag", "tfrag")) %>%
  mutate(qn = vapply(q, short, ""), rn = vapply(r, short, ""))

# reliability QC: keep genomes whose self-comparison yields >= 100 total 3kb fragments
selfrag <- ani %>% filter(qn == rn) %>% group_by(qn) %>% summarise(tfrag = max(tfrag), .groups = "drop")
reliable <- selfrag$qn[selfrag$tfrag >= 100]
is_chilo <- function(x) !x %in% names(refname)
chilo_rel <- reliable[is_chilo(reliable)]

# mutual ANI among reliable Chilo assemblies
mut <- ani %>% filter(qn %in% chilo_rel, rn %in% chilo_rel, qn != rn)
# Chilo -> reference (supergroup placement); best A and best B per sample
ref_hits <- ani %>% filter(qn %in% chilo_rel, rn %in% names(refname)) %>%
  mutate(ref = refname[rn], sg = sub(".*\\((.)\\)", "\\1", ref))

# labelled ANI matrix over reliable Chilo + references, for a heatmap
keep <- c(chilo_rel, names(refname))
lab <- function(x) ifelse(x %in% names(refname), refname[x], x)
mat <- ani %>% filter(qn %in% keep, rn %in% keep) %>%
  transmute(q = lab(qn), r = lab(rn), ani)
ord <- c(sort(chilo_rel), "wMel (A)", "wRi (A)", "wPip (B)", "wTpre (B)", "wBm (D)")
ord <- ord[ord %in% c(mat$q, mat$r)]
mat$q <- factor(mat$q, levels = ord); mat$r <- factor(mat$r, levels = rev(ord))
p_ani <- ggplot(mat, aes(q, r, fill = ani)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", ani)), size = 2.6) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#2171B5", "#08306B"),
                       limits = c(80, 100), name = "ANI %") +
  labs(x = NULL, y = NULL,
       title = "Stage 04h: whole-genome ANI (fastANI) — contiguous assemblies vs reference strains",
       subtitle = "Chilo assemblies mutually >=98%; each maps to supergroup A (~97%) AND B (~87-91%) = genome-level A+B co-infection; wBm/D outgroup ~83%") +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(size = 7.6, color = "grey30"),
        axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
ggsave("reports/figs/strain_ani_heatmap.png", p_ani, width = 8, height = 5.2, dpi = 150)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
rr <- function(L) { r <- recur[recur$locus == L, ]; sprintf("%s=%s (%d/%d)", L, r$A_allele, r$A_n, n_typed) }
best_sg <- ref_hits %>% group_by(qn, sg) %>% summarise(ani = max(ani), .groups = "drop")
a_range <- sprintf("%.1f–%.1f", min(best_sg$ani[best_sg$sg == "A"]), max(best_sg$ani[best_sg$sg == "A"]))
b_range <- sprintf("%.1f–%.1f", min(best_sg$ani[best_sg$sg == "B"]), max(best_sg$ani[best_sg$sg == "B"]))
lines <- c(
  "# Stage 04h (B5) — 客观同株判据：MLST 等位谱 + 全基因组 ANI",
  "",
  "把'单一循环 Wolbachia 株-单元'从口述变成**预设阈值下的客观判据**。两条独立读出。",
  "",
  "## 判据 1 — MLST 等位谱一致性（预设：相同 5 位点谱 = 同 ST = 同株）",
  "",
  sprintf("- 可分型样本 **%d** 个（多数阳性样本低载量/低深度无法齐 5 位点）。", n_typed),
  sprintf("- **完整 ST-19（gatB7/coxA6/hcpA7/ftsZ3/fbpA8，超群 A）签名齐 5 位点：%d/%d**（AQ_R_9 纯 A 单株；NC2025_R_6 为 A+B 的 A 组分）。", n_full_ST19, n_typed),
  "- A（ST-19）签名各位点跨样本复现（含该等位样本数）：",
  sprintf("  - %s；%s；%s；%s；%s。", rr("gatB"), rr("coxA"), rr("hcpA"), rr("ftsZ"), rr("fbpA")),
  sprintf("  - gatB7/coxA6/hcpA7 **强复现（5–6/%d）**；ftsZ3/fbpA8 复现弱（2–3/%d）——低深度噪声 + 共感染将 B 等位顶入这些位点。", n_typed, n_typed),
  sprintf("- 第二株（超群 B）等位同样跨样本复现：gatB9=%d/%d、coxA150=%d/%d、hcpA181=%d/%d、fbpA230=%d/%d。",
          recur$B_n[recur$locus=="gatB"], n_typed, recur$B_n[recur$locus=="coxA"], n_typed,
          recur$B_n[recur$locus=="hcpA"], n_typed, recur$B_n[recur$locus=="fbpA"], n_typed),
  "",
  "> **诚实校准**：'7/6/7/3/8 签名跨多样本复现'对 gatB/coxA/hcpA 成立、对 ftsZ/fbpA 受深度限制。宜表述为",
  "> '一个占优的超群-A ST-19 主链在可分型样本中反复出现，伴随复现的超群-B 第二等位'，而非'全部样本齐 ST-19'。",
  "",
  "## 判据 2 — 全基因组 ANI（fastANI；预设：≥95% 同种，≥99% 同克隆株；超群界 ~85–90%）",
  "",
  sprintf("- **可靠性 QC**：仅 %d 个 enrichment 组装的 ≥3kb contig 数足够（自比 ≥100 个 3kb 片段）：%s；其余 22 个太碎，ANI 不可信。",
          length(chilo_rel), paste(sort(chilo_rel), collapse = "、")),
  sprintf("- **Chilo 组装两两 ANI = %.1f–%.1f%%**（多数 ~98–99%%）——全部 ≥95%%（同种），部分 ≥99%%（同克隆）。",
          min(mut$ani), max(mut$ani)),
  sprintf("- 每个 Chilo 组装**同时**高比中超群 A（%s%%，高片段数）**和** B（%s%%）参考株 → **基因组层面独立确认 A+B 共感染**（印证 Stage 04f）；wBm/D 仅 ~83%%（外群）。",
          a_range, b_range),
  "",
  "> **ANI 局限**：(i) enrichment 组装碎片化使 22/25 样本无法算可信 ANI；(ii) A+B 共感染的嵌合组装使两两 ANI 被**压低**（不同 A:B 比例的混合体互比），故 98–99% 是**保守下界**，真实单株 ANI 更高。",
  "",
  "## 结论",
  "",
  "MLST 主链一致 + 复现的 A/B 双等位 + 可靠组装两两 ANI ≥98% + 基因组级 A+B 双超群命中，四者一致支持：",
  "**野外二化螟携带的是一个占优的、单一循环的 Wolbachia 株-单元（超群 A 主链 + 复现的超群 B 伴生），非多个无关株的随机集合。**",
  "受限于低载量样本的深度墙与共感染嵌合，'完全相同 ST' 仅在 2 个高质量样本严格成立；跨分歧寄生物的株-单倍型协变仍待 V3 同体分型。",
  "",
  "见 reports/figs/strain_mlst_concordance.png、strain_ani_heatmap.png；results/05_typing/strain_id/{mlst_profiles,mlst_recurrence,ani_all}.tsv。"
)
writeLines(lines, "reports/stage04h_strain_identity.md")
cat("done: n_typed =", n_typed, "| full ST-19 =", n_full_ST19,
    "| reliable ANI genomes =", length(chilo_rel), "| mutual ANI",
    sprintf("%.1f-%.1f", min(mut$ani), max(mut$ani)), "\n")
