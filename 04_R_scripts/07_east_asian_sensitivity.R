# ============================================================
# 07_east_asian_sensitivity.R
# Project: 01_MR_GastricCancer
#
# Prospective analysis after OSF registration
# East Asian cross-ancestry sensitivity outcome:
# GCST90018629
# ============================================================

library(TwoSampleMR)

source("04_R_scripts/00_functions.R")


# ============================================================
# 1. Output folders
# ============================================================

dir.create(
  "02_intermediate_data/02_outcome/EastAsian",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "05_results/04_east_asian_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "06_figures/04_east_asian_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Load PRIMARY microbiome instruments
#    P <= 5e-6
# ============================================================

veillonella_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

veillonellaceae_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_primary_p5e6.rds"
)

nrow(veillonella_iv)
nrow(veillonellaceae_iv)
# ============================================================
# 3. Prespecified SNP list
# ============================================================

target_snps_eas <- unique(
  c(
    veillonella_iv$SNP,
    veillonellaceae_iv$SNP
  )
)

length(target_snps_eas)
target_snps_eas
# ============================================================
# 4. East Asian gastric cancer VCF
# ============================================================

file_eas <-
  "01_raw_data/02_outcome/ebi-a-GCST90018629.vcf.gz"

file.exists(file_eas)
# ============================================================
# 5. Extract prespecified IVs from GCST90018629
# ============================================================

eas_result <- extract_outcome_snps(
  vcf_file = file_eas,
  target_snps = target_snps_eas,
  chunk_size = 50000
)

eas_raw <- eas_result$data
# ============================================================
# 6. SNP matching QC
# ============================================================

eas_result$n_requested
eas_result$n_found
eas_result$n_missing
eas_result$missing
# ============================================================
# 7. East Asian outcome data QC
# ============================================================

nrow(eas_raw)

sum(is.na(eas_raw$beta))
sum(is.na(eas_raw$se))
sum(is.na(eas_raw$pval))
sum(is.na(eas_raw$eaf))

sum(duplicated(eas_raw$SNP))

eas_raw
# 8.保存缺失SNP/missing snp
write.csv(
  data.frame(SNP = eas_result$missing),
  "02_intermediate_data/02_outcome/EastAsian/GCST90018629_missing_prespecified_IVs.csv",
  row.names = FALSE
)
#保存剩余snp
saveRDS(
  eas_raw,
  "02_intermediate_data/02_outcome/EastAsian/GCST90018629_extracted_IVs.rds"
)

write.csv(
  eas_raw,
  "02_intermediate_data/02_outcome/EastAsian/GCST90018629_extracted_IVs.csv",
  row.names = FALSE
)
# ============================================================
# 9. Format East Asian outcome for TwoSampleMR
# ============================================================

eas_format <- eas_raw

eas_format$Phenotype <- "East Asian gastric cancer"

eas_outcome <- TwoSampleMR::format_data(
  eas_format,
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

nrow(eas_outcome)
# ============================================================
# 10. Separate outcome SNPs by exposure
# ============================================================

eas_veillonella <- eas_outcome[
  eas_outcome$SNP %in% veillonella_iv$SNP,
]

eas_veillonellaceae <- eas_outcome[
  eas_outcome$SNP %in% veillonellaceae_iv$SNP,
]

nrow(eas_veillonella)
nrow(eas_veillonellaceae)
# ============================================================
# 11. Harmonisation
# action = 3: conservatively remove palindromic SNPs
# ============================================================

harm_veillonella_eas <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonella_iv,
  outcome_dat = eas_veillonella,
  action = 3
)

harm_veillonellaceae_eas <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonellaceae_iv,
  outcome_dat = eas_veillonellaceae,
  action = 3
)
# ============================================================
# 12. Harmonisation QC
# ============================================================

table(harm_veillonella_eas$mr_keep)
table(harm_veillonellaceae_eas$mr_keep)

harm_veillonella_eas[
  harm_veillonella_eas$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]

harm_veillonellaceae_eas[
  harm_veillonellaceae_eas$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]
#生成MR-ready数据
veillonella_eas_ready <- subset(
  harm_veillonella_eas,
  mr_keep == TRUE
)

veillonellaceae_eas_ready <- subset(
  harm_veillonellaceae_eas,
  mr_keep == TRUE
)

nrow(veillonella_eas_ready)
nrow(veillonellaceae_eas_ready)

sum(is.na(veillonella_eas_ready$beta.exposure))
sum(is.na(veillonella_eas_ready$beta.outcome))

sum(is.na(veillonellaceae_eas_ready$beta.exposure))
sum(is.na(veillonellaceae_eas_ready$beta.outcome))
# ============================================================
# 13. MR methods
# ============================================================

mr_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)


# ============================================================
# 14. East Asian MR: Veillonella -> gastric cancer
# ============================================================

mr_veillonella_eas <- TwoSampleMR::mr(
  veillonella_eas_ready,
  method_list = mr_methods
)

mr_veillonella_eas


# ============================================================
# 15. East Asian MR: Veillonellaceae -> gastric cancer
# ============================================================

mr_veillonellaceae_eas <- TwoSampleMR::mr(
  veillonellaceae_eas_ready,
  method_list = mr_methods
)

mr_veillonellaceae_eas
#转换OR值
# ============================================================
# 16. Convert to odds ratios
# ============================================================

mr_veillonella_eas_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonella_eas
  )

mr_veillonellaceae_eas_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonellaceae_eas
  )

mr_veillonella_eas_or
mr_veillonellaceae_eas_or
# 17 Heterogeneity
het_veillonella_eas <- TwoSampleMR::mr_heterogeneity(
  veillonella_eas_ready
)

het_veillonellaceae_eas <- TwoSampleMR::mr_heterogeneity(
  veillonellaceae_eas_ready
)

het_veillonella_eas
het_veillonellaceae_eas


# 18 MR-Egger intercept
pleio_veillonella_eas <- TwoSampleMR::mr_pleiotropy_test(
  veillonella_eas_ready
)

pleio_veillonellaceae_eas <- TwoSampleMR::mr_pleiotropy_test(
  veillonellaceae_eas_ready
)

pleio_veillonella_eas
pleio_veillonellaceae_eas
# ============================================================
# 19. Leave-one-out
# ============================================================

loo_veillonella_eas <- TwoSampleMR::mr_leaveoneout(
  veillonella_eas_ready
)

loo_veillonellaceae_eas <- TwoSampleMR::mr_leaveoneout(
  veillonellaceae_eas_ready
)

loo_veillonella_eas
loo_veillonellaceae_eas


# ============================================================
# 20. Single-SNP analysis
# ============================================================

single_veillonella_eas <- TwoSampleMR::mr_singlesnp(
  veillonella_eas_ready
)

single_veillonellaceae_eas <- TwoSampleMR::mr_singlesnp(
  veillonellaceae_eas_ready
)

single_veillonella_eas
single_veillonellaceae_eas
# ============================================================
# 21. MR-PRESSO
# ============================================================

library(MRPRESSO)

presso_veillonella_eas <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonella_eas_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonellaceae_eas <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonellaceae_eas_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonella_eas
presso_veillonellaceae_eas
# ============================================================
# 22. Save East Asian MR results
# ============================================================

write.csv(
  mr_veillonella_eas_or,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_MR.csv",
  row.names = FALSE
)

write.csv(
  mr_veillonellaceae_eas_or,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_MR.csv",
  row.names = FALSE
)

write.csv(
  het_veillonella_eas,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  het_veillonellaceae_eas,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonella_eas,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonellaceae_eas,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonella_eas,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonellaceae_eas,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  single_veillonella_eas,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  single_veillonellaceae_eas,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_singleSNP.csv",
  row.names = FALSE
)

saveRDS(
  presso_veillonella_eas,
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_MRPRESSO.rds"
)

saveRDS(
  presso_veillonellaceae_eas,
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_MRPRESSO.rds"
)

capture.output(
  presso_veillonella_eas,
  file = "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_MRPRESSO.txt"
)

capture.output(
  presso_veillonellaceae_eas,
  file = "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_MRPRESSO.txt"
)
# ============================================================
# 23. Save East Asian figures
# ============================================================

forest_veillonella_eas <-
  TwoSampleMR::mr_forest_plot(single_veillonella_eas)

forest_veillonellaceae_eas <-
  TwoSampleMR::mr_forest_plot(single_veillonellaceae_eas)

loo_plot_veillonella_eas <-
  TwoSampleMR::mr_leaveoneout_plot(loo_veillonella_eas)

loo_plot_veillonellaceae_eas <-
  TwoSampleMR::mr_leaveoneout_plot(loo_veillonellaceae_eas)

ggplot2::ggsave(
  "06_figures/04_east_asian_sensitivity/Veillonella_GCST90018629_singleSNP_forest.pdf",
  plot = forest_veillonella_eas[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_singleSNP_forest.pdf",
  plot = forest_veillonellaceae_eas[[1]],
  width = 7,
  height = 7
)

ggplot2::ggsave(
  "06_figures/04_east_asian_sensitivity/Veillonella_GCST90018629_leaveoneout.pdf",
  plot = loo_plot_veillonella_eas[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_leaveoneout.pdf",
  plot = loo_plot_veillonellaceae_eas[[1]],
  width = 7,
  height = 7
)
#检查
list.files("05_results/04_east_asian_sensitivity")
list.files("06_figures/04_east_asian_sensitivity")
