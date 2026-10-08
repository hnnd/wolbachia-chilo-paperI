# Analysis code and figure source data — Paper I

**Manuscript.** *Wolbachia* detected in field populations of the rice stem borer *Chilo
suppressalis* originates from its endoparasitoids: signal attribution, lineage specificity,
and two co-occurring genomes (referred to here as Paper I).

**Version.** v1.0.0 — 2026-10-08

**Authors.** Shiqing Yuan, Pinbo Liu, Meng Wang, Jia Li, Lijuan Yin, Yunsheng Wang,
Yongmo Wang.

**Contact.** Yunsheng Wang (wangys@hunau.edu.cn), Yongmo Wang (ymwang@mail.hzau.edu.cn).

**License.** Creative Commons Attribution 4.0 International (CC BY 4.0), see `LICENSE`.

---

## Contents

| Path | Content |
|---|---|
| `code/workflow/` | The analysis pipeline: numbered SLURM stage scripts (`NN_*.sh`), per-sample array workers (`jobs/`), shared shell library (`lib/`), figure scripts (`figures/`), competitive-mapping scripts (`mapping/`), primer design (`primers/`) |
| `code/config/` | `config.yaml` (all thresholds and paths), `samples.tsv` (the 247-library manifest), `build_samples.sh` (rebuilds `samples.tsv`) |
| `source_data/` | Figure source tables and key diagnostic tables analysed for the paper (`data/tables/` of the working repository) |
| `figures/` | The final main figures (`Figure_1.pdf` … `Figure_5.pdf`) and Supporting Information figures (`Figure_S1.pdf` … `Figure_S10.pdf`), exactly as submitted |
| `docs/` | Design document (`PAPER_I_DESIGN.md`), the v5 assembly audit (`2026-09-30-v5-audit-and-deposit-fix.md`), and the two method reports (`LANE_RETROSPECTIVE_REPORT.md`, `HOST_PLANT_PARASITOID_ANALYSIS.md`) |
| `MANIFEST.sha256` | SHA-256 checksums of every file in this archive |

## Study data (public repositories)

- **210 newly sequenced libraries** — NCBI BioProject **PRJNA1527622**
  (BioSample SAMN63091504–SAMN63091713; one per library, listed in
  `source_data/TableS1_per_library_criteria.csv`)
- **37 public whole-body datasets** — BioProject **PRJNA1011001**
- **PacBio HiFi source for the two genomes** — BioProject **PRJNA1480045**
  (source BioSample SAMN61026177)
- **Wolbachia genomes** — *w*CchiA: SAMN63092028, GenBank WGS **JCEMCJ000000000**;
  *w*CchiB: SAMN63092029, GenBank WGS **JCEMCI000000000** (submission SUB16507602)
- ***wsp* Sanger consensus sequences** — GenBank **QB079202** (JSR-1), **QB079203** (JSR-2)
  (submission SUB16514567)

## Pipeline overview

The analysis is organized as numbered stages, each a thin submitter plus a SLURM array
worker, with per-stage `.done` sentinels so any stage can be re-run incrementally. All
thresholds live in `code/config/config.yaml`; the random seed is fixed at 42.

| Stage | Purpose | Output |
|---|---|---|
| 00a/00b | conda environments; host, *Wolbachia* (4 strains), parasitoid and Kraken2 references | `refs/` |
| 01 | read QC, host dehosting (`minibwa map -x sr`) | `results/02_nonhost/` |
| 02 / 02b | per-library profiling vs host, each *Wolbachia* strain (highest-breadth kept), parasitoid panel; positivity calls | `results/03_profile/per_sample_profile.tsv` |
| 02c | Kraken2 per-sample taxonomy (sidecar) | `results/03_profile/*.kraken.report` |
| 03 | within-sample co-occurrence analysis | `reports/stage03*.md` |
| 03b–03h | dose–response, mito dual-source, census, de novo parasitoid ID, threshold sensitivity, multivariable (Firth) model | `data/tables/`, `reports/stage03*.md` |
| 04 / 04b | SPAdes metagenome assembly of positives; MLST/*wsp* typing; CheckM2 | `results/04_assembly/`, `results/05_typing/` |
| 04c–04h | phylogenomics, *wsp* tree, consensus typing, co-infection, co-phylogeny, strain identity | `reports/stage04*.md` |
| 17 | HiFi read coverage of the deposited genomes | `results/17_ncbi_cov/coverage_summary.tsv` |

Re-running requires a SLURM cluster and the external references (`code/workflow/lib/fetch_*.sh`
download the public ones). Paths in `code/config/config.yaml` and `code/config/samples.tsv`
are absolute paths from the original cluster and must be adapted; the FASTQ files themselves
are in the BioProjects above, not in this archive.

## Software

fastp 1.3.5 · minibwa 0.3 (r391, commit a679d17) · Kraken2 2.1.3 · SAMtools 1.21 · seqkit
2.13.0 · csvtk 0.37.0 · SPAdes/metaSPAdes 4.0.0 · CheckM2 1.1.0 (DIAMOND 2.1.11) · mlst
2.35.0 · MAFFT 7.526 · IQ-TREE 3.1.2 · fastANI 1.34 · BLAST+ 2.16.0 · MEGAHIT 1.2.9 ·
RagTag 2.1.0 · Prodigal 2.6.3 · barrnap 1.10.5 · tRNAscan-SE 2.0.13 · minimap2 2.28-r1209
and 2.31-r1302 · Biopython 1.87 · R 4.5.3 (logistf 1.26.1, ggplot2 4.0.3, dplyr 1.2.1) ·
Python 3.12.3 (pandas, matplotlib). Full version list: §2.12 of the manuscript.

## Figure map

Final figures are in `figures/`; source data are in `source_data/`. Some figure-composition
scripts were run on the analysis cluster or in an analysis session and are **not** part of
this archive; where that is the case it is stated below. The manuscript's figure legends and
`docs/PAPER_I_DESIGN.md` document what each panel shows.

| Figure | Composition script in this archive | Source data in `source_data/` |
|---|---|---|
| Figure 1 | — (not included) | `figI1_source_data.csv`, `figI1_dose_bins.csv`, `TableS1_per_library_criteria.csv` |
| Figure 2 | `code/workflow/figures/make_figI2.py` | `figI2_ani_matrix.csv`, `figI2_genome_quality.csv`, `figI2_tree_distances.csv`, `figI2_btype_identities.csv`, `figI2_wsp_*` |
| Figure 3 | — (not included) | `figI3_taxon_summary.csv`, `figI3_lineage_carriage.csv`, `figI3_site_composition.csv`, `cotesia_verification_status.csv`, `denovo_lowload_four.csv` |
| Figure 4 | — (not included) | `host_plant_summary.csv`, `parasitoid_taxon_by_host_plant.csv`, `cotesia_by_site.csv`, `per_site_paired.csv`; method report in `docs/` |
| Figure 5 | — (composite of microscopy images; not included) | `figI5_contrast_check.csv`, `figI5_image_adjustment.csv` |
| Figure S1 | `code/workflow/figures/build_figI_S1.py` | `TableS1_per_library_criteria.csv` |
| Figure S2 | `code/workflow/figures/build_figI_S2.py` | `coverage_uniformity*.csv`, `wolb_diag_*.csv` |
| Figure S3 | — (composition script not included; the threshold-scan analysis is implemented in `code/workflow/03g_threshold_sensitivity.R`) | threshold scan derived from `TableS1_per_library_criteria.csv` |
| Figure S4 | — (composition script not included; the multivariable model is implemented in `code/workflow/03h_multivariate.R`; the submitted panel reports the 33-case refit) | `multivariate_coef_33.tsv` |
| Figure S5 | — (not included) | see `docs/LANE_RETROSPECTIVE_REPORT.md` |
| Figure S6 | `code/workflow/figures/build_figI_S6.py` | `figI2_genome_quality.csv`, `figI2_removed_protein_identity.csv`, `unplaced_protein_homology.csv` |
| Figure S7 | `code/workflow/figures/build_figI_S7.py` | counts documented inline in the script |
| Figure S8 | — (not included) | `negative_sample_decomposition.csv`, `coverage_uniformity*.csv`, `wolb_diag_calibration.csv` |
| Figure S9 | data: `code/workflow/mapping/ab_competitive_v4.sh`; composition script not included | `ab_competitive_v4.tsv`, `TableS3_strain_assignment.csv`, `bin_cocoverage_matrix.csv` |
| Figure S10 | — (not included) | `figI4_window_identity.csv`, `figI4_window_significance.csv`, `figI4_recomb_pergene.csv`, `figI4_recomb_candidates.csv`, `binB_window_identity.csv` |

Reproducibility notes:

- Earlier-iteration scripts are included for provenance only and cannot reproduce the
  submitted panels: `code/workflow/figures/build_figs.py` asserts the earlier 31-positive
  call set (its output set predates the final reclassification), as do the
  `fig1_cooccurrence.R`–`fig4_evidenceB.R` scripts. The figure map above routes to the
  scripts written for the submitted version.
- `make_figI2.py` reads its inputs (core-gene tree, fastANI table, whole-genome alignment)
  from analysis-session paths recorded in the script; those inputs are regenerated by the
  cluster stages described in §2.5 of the manuscript and the design document.
- Several composition scripts under `code/workflow/figures/` were written against the
  working repository layout (`results/`, `reports/figs/`) and expect those directories to
  exist; they document the exact plotting choices and can be re-pointed at the files in
  `source_data/`.
- `samples.tsv` carries the original cluster FASTQ paths; after downloading from the
  BioProjects, run `build_samples.sh` (or edit the paths) to point at the local copies.

## Citation

If you use this code or these tables, please cite the manuscript and this archive:

> Yuan, S., Liu, P., Wang, M., Li, J., Yin, L., Wang, Y., & Wang, Y. (2026).
> Analysis code and figure source data for "Wolbachia detected in field populations of the
> rice stem borer Chilo suppressalis originates from its endoparasitoids" (v1.0.0)
> [Computer software]. Zenodo. https://doi.org/—— (DOI assigned on publication)

A machine-readable record is provided in `CITATION.cff`.
