# ============================================================
# 08_reverse_MR.R
# Project: 01_MR_GastricCancer
#
# Prospective reverse-direction MR after OSF registration
# Exposure: genetic liability to gastric cancer
# Outcomes:
#   1. Veillonella
#   2. Veillonellaceae
# ============================================================

library(TwoSampleMR) 

source("04_R_scripts/00_functions.R")


# ============================================================
# 1. Output folders
# ============================================================

dir.create(
  "02_intermediate_data/03_reverse_MR",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "05_results/05_reverse_MR",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "06_figures/05_reverse_MR",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Check OpenGWAS token
# ============================================================

nzchar(Sys.getenv("OPENGWAS_JWT"))
# ============================================================
# 3. Gastric cancer instruments
#    P < 5e-8
#    LD r2 < 0.001
#    clumping window = 10,000 kb
# ============================================================

gc_reverse_iv <- TwoSampleMR::extract_instruments(
  outcomes = "ebi-a-GCST90018849",
  p1 = 5e-8,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)
nrow(gc_reverse_iv)

gc_reverse_iv[, c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "effect_allele.exposure",
  "other_allele.exposure"
)]
# ============================================================
# 4. Instrument strength
# ============================================================

gc_reverse_iv$F_stat <-
  (gc_reverse_iv$beta.exposure /
     gc_reverse_iv$se.exposure)^2

gc_reverse_iv[, c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "F_stat"
)]

summary(gc_reverse_iv$F_stat)
#按预注册规则排除
gc_reverse_iv <- subset(
  gc_reverse_iv,
  F_stat > 10
)

nrow(gc_reverse_iv)
gc_reverse_iv$SNP
#保存
saveRDS(
  gc_reverse_iv,
  "02_intermediate_data/03_reverse_MR/GCST90018849_reverse_IVs.rds"
)

write.csv(
  gc_reverse_iv,
  "02_intermediate_data/03_reverse_MR/GCST90018849_reverse_IVs.csv",
  row.names = FALSE
)
# ============================================================
# 5. Define reverse-MR microbiome outcome files
# ============================================================

file_veillonella_rev <-
  "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90017088.vcf.gz"

file_veillonellaceae_rev <-
  "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90016956.vcf.gz"

file.exists(file_veillonella_rev)
file.exists(file_veillonellaceae_rev)
# ============================================================
# 6. Extract GC instruments from Veillonella outcome GWAS
# ============================================================

rev_veillonella_result <- extract_outcome_snps(
  vcf_file = file_veillonella_rev,
  target_snps = gc_reverse_iv$SNP,
  chunk_size = 50000
)

rev_veillonella_raw <- rev_veillonella_result$data


rev_veillonella_result$n_requested
rev_veillonella_result$n_found
rev_veillonella_result$n_missing
rev_veillonella_result$missing
# ============================================================
# 7. Extract GC instruments from Veillonellaceae outcome GWAS
# ============================================================

rev_veillonellaceae_result <- extract_outcome_snps(
  vcf_file = file_veillonellaceae_rev,
  target_snps = gc_reverse_iv$SNP,
  chunk_size = 50000
)

rev_veillonellaceae_raw <- rev_veillonellaceae_result$data


rev_veillonellaceae_result$n_requested
rev_veillonellaceae_result$n_found
rev_veillonellaceae_result$n_missing
rev_veillonellaceae_result$missing
# ============================================================
# 8. Reverse outcome QC
# ============================================================

# Veillonella
nrow(rev_veillonella_raw)

sum(is.na(rev_veillonella_raw$beta))
sum(is.na(rev_veillonella_raw$se))
sum(is.na(rev_veillonella_raw$pval))
sum(is.na(rev_veillonella_raw$eaf))
sum(duplicated(rev_veillonella_raw$SNP))


# Veillonellaceae
nrow(rev_veillonellaceae_raw)

sum(is.na(rev_veillonellaceae_raw$beta))
sum(is.na(rev_veillonellaceae_raw$se))
sum(is.na(rev_veillonellaceae_raw$pval))
sum(is.na(rev_veillonellaceae_raw$eaf))
sum(duplicated(rev_veillonellaceae_raw$SNP))
# ============================================================
# 9. Save unavailable reverse-MR instruments
# ============================================================

write.csv(
  data.frame(SNP = rev_veillonella_result$missing),
  "02_intermediate_data/03_reverse_MR/Veillonella_missing_GC_IVs.csv",
  row.names = FALSE
)

write.csv(
  data.frame(SNP = rev_veillonellaceae_result$missing),
  "02_intermediate_data/03_reverse_MR/Veillonellaceae_missing_GC_IVs.csv",
  row.names = FALSE
)
saveRDS(
  rev_veillonella_raw,
  "02_intermediate_data/03_reverse_MR/Veillonella_reverse_outcome_raw.rds"
)

saveRDS(
  rev_veillonellaceae_raw,
  "02_intermediate_data/03_reverse_MR/Veillonellaceae_reverse_outcome_raw.rds"
)
# ============================================================
# 10. Format microbiome outcomes
# ============================================================

rev_veillonella_format <- rev_veillonella_raw
rev_veillonella_format$Phenotype <- "Veillonella"

rev_veillonella_outcome <- TwoSampleMR::format_data(
  rev_veillonella_format,
  type = "outcome",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)


rev_veillonellaceae_format <- rev_veillonellaceae_raw
rev_veillonellaceae_format$Phenotype <- "Veillonellaceae"

rev_veillonellaceae_outcome <- TwoSampleMR::format_data(
  rev_veillonellaceae_format,
  type = "outcome",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

nrow(rev_veillonella_outcome)
nrow(rev_veillonellaceae_outcome)
# ============================================================
# 11. Reverse harmonisation
# ============================================================

harm_gc_veillonella <- TwoSampleMR::harmonise_data(
  exposure_dat = gc_reverse_iv,
  outcome_dat = rev_veillonella_outcome,
  action = 3
)

harm_gc_veillonellaceae <- TwoSampleMR::harmonise_data(
  exposure_dat = gc_reverse_iv,
  outcome_dat = rev_veillonellaceae_outcome,
  action = 3
)
# ============================================================
# 12. Harmonisation QC
# ============================================================

table(harm_gc_veillonella$mr_keep)
table(harm_gc_veillonellaceae$mr_keep)

harm_gc_veillonella[
  harm_gc_veillonella$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]

harm_gc_veillonellaceae[
  harm_gc_veillonellaceae$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]

gc_veillonella_ready <- subset(
  harm_gc_veillonella,
  mr_keep == TRUE
)

gc_veillonellaceae_ready <- subset(
  harm_gc_veillonellaceae,
  mr_keep == TRUE
)

nrow(gc_veillonella_ready)
nrow(gc_veillonellaceae_ready)

sum(is.na(gc_veillonella_ready$beta.exposure))
sum(is.na(gc_veillonella_ready$beta.outcome))

sum(is.na(gc_veillonellaceae_ready$beta.exposure))
sum(is.na(gc_veillonellaceae_ready$beta.outcome))
# ============================================================
# 13. Reverse MR methods
# ============================================================

mr_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)


# ============================================================
# 14. Gastric cancer -> Veillonella
# ============================================================

mr_gc_veillonella <- TwoSampleMR::mr(
  gc_veillonella_ready,
  method_list = mr_methods
)

mr_gc_veillonella


# ============================================================
# 15. Gastric cancer -> Veillonellaceae
# ============================================================

mr_gc_veillonellaceae <- TwoSampleMR::mr(
  gc_veillonellaceae_ready,
  method_list = mr_methods
)

mr_gc_veillonellaceae
# ============================================================
# 16. Reverse MR heterogeneity
# ============================================================

het_gc_veillonella <- TwoSampleMR::mr_heterogeneity(
  gc_veillonella_ready
)

het_gc_veillonellaceae <- TwoSampleMR::mr_heterogeneity(
  gc_veillonellaceae_ready
)

het_gc_veillonella
het_gc_veillonellaceae


# ============================================================
# 17. MR-Egger intercept
# ============================================================

pleio_gc_veillonella <- TwoSampleMR::mr_pleiotropy_test(
  gc_veillonella_ready
)

pleio_gc_veillonellaceae <- TwoSampleMR::mr_pleiotropy_test(
  gc_veillonellaceae_ready
)

pleio_gc_veillonella
pleio_gc_veillonellaceae
# ============================================================
# 18. Leave-one-out
# ============================================================

loo_gc_veillonella <- TwoSampleMR::mr_leaveoneout(
  gc_veillonella_ready
)

loo_gc_veillonellaceae <- TwoSampleMR::mr_leaveoneout(
  gc_veillonellaceae_ready
)


# ============================================================
# 19. Single-SNP
# ============================================================

single_gc_veillonella <- TwoSampleMR::mr_singlesnp(
  gc_veillonella_ready
)

single_gc_veillonellaceae <- TwoSampleMR::mr_singlesnp(
  gc_veillonellaceae_ready
)# ============================================================
# 18. Leave-one-out
# ============================================================

loo_gc_veillonella <- TwoSampleMR::mr_leaveoneout(
  gc_veillonella_ready
)

loo_gc_veillonellaceae <- TwoSampleMR::mr_leaveoneout(
  gc_veillonellaceae_ready
)


# ============================================================
# 19. Single-SNP
# ============================================================

single_gc_veillonella <- TwoSampleMR::mr_singlesnp(
  gc_veillonella_ready
)

single_gc_veillonellaceae <- TwoSampleMR::mr_singlesnp(
  gc_veillonellaceae_ready
)
# ============================================================
# 20. MR-PRESSO
# ============================================================

library(MRPRESSO)

presso_gc_veillonella <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = gc_veillonella_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_gc_veillonellaceae <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = gc_veillonellaceae_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_gc_veillonella
presso_gc_veillonellaceae
# ============================================================
# 21. Save reverse MR results
# ============================================================

write.csv(
  mr_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_MR.csv",
  row.names = FALSE
)

write.csv(
  mr_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_MR.csv",
  row.names = FALSE
)

write.csv(
  het_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  het_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  pleio_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  pleio_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  loo_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  loo_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  single_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  single_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_singleSNP.csv",
  row.names = FALSE
)

saveRDS(
  presso_gc_veillonella,
  "05_results/05_reverse_MR/GC_to_Veillonella_MRPRESSO.rds"
)

saveRDS(
  presso_gc_veillonellaceae,
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_MRPRESSO.rds"
)

capture.output(
  presso_gc_veillonella,
  file = "05_results/05_reverse_MR/GC_to_Veillonella_MRPRESSO.txt"
)

capture.output(
  presso_gc_veillonellaceae,
  file = "05_results/05_reverse_MR/GC_to_Veillonellaceae_MRPRESSO.txt"
)
# ============================================================
# 22. Reverse MR figures
# ============================================================

forest_gc_veillonella <-
  TwoSampleMR::mr_forest_plot(single_gc_veillonella)

forest_gc_veillonellaceae <-
  TwoSampleMR::mr_forest_plot(single_gc_veillonellaceae)

loo_plot_gc_veillonella <-
  TwoSampleMR::mr_leaveoneout_plot(loo_gc_veillonella)

loo_plot_gc_veillonellaceae <-
  TwoSampleMR::mr_leaveoneout_plot(loo_gc_veillonellaceae)

ggplot2::ggsave(
  "06_figures/05_reverse_MR/GC_to_Veillonella_singleSNP_forest.pdf",
  plot = forest_gc_veillonella[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/05_reverse_MR/GC_to_Veillonellaceae_singleSNP_forest.pdf",
  plot = forest_gc_veillonellaceae[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/05_reverse_MR/GC_to_Veillonella_leaveoneout.pdf",
  plot = loo_plot_gc_veillonella[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/05_reverse_MR/GC_to_Veillonellaceae_leaveoneout.pdf",
  plot = loo_plot_gc_veillonellaceae[[1]],
  width = 7,
  height = 6
)
# ============================================================
# Recreate reverse MR plotting objects
# ============================================================

forest_gc_veillonella <-
  TwoSampleMR::mr_forest_plot(single_gc_veillonella)

forest_gc_veillonellaceae <-
  TwoSampleMR::mr_forest_plot(single_gc_veillonellaceae)

loo_plot_gc_veillonella <-
  TwoSampleMR::mr_leaveoneout_plot(loo_gc_veillonella)

loo_plot_gc_veillonellaceae <-
  TwoSampleMR::mr_leaveoneout_plot(loo_gc_veillonellaceae)
#检查
list.files("06_figures/05_reverse_MR")
list.files("05_results/05_reverse_MR")
