# ============================================================
# 06_finngen_sensitivity.R
# Project: 01_MR_GastricCancer
#
# Prospective analysis after OSF registration
# Secondary European sensitivity outcome:
# FinnGen gastric cancer
# ============================================================

library(TwoSampleMR)

source("04_R_scripts/00_functions.R")


# ============================================================
# 1. Output folders
# ============================================================

dir.create(
  "02_intermediate_data/02_outcome/FinnGen",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "05_results/03_finngen_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "06_figures/03_finngen_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Load PRIMARY microbiome IVs
#    Instrument threshold: P <= 5e-6
# ============================================================

veillonella_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

veillonellaceae_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_primary_p5e6.rds"
)

# Check
nrow(veillonella_iv)
nrow(veillonellaceae_iv)

veillonella_iv$SNP
veillonellaceae_iv$SNP


# ============================================================
# 3. Prespecified SNP list
# ============================================================

target_snps_finn <- unique(
  c(
    veillonella_iv$SNP,
    veillonellaceae_iv$SNP
  )
)

length(target_snps_finn)
target_snps_finn
#检查FinnGen文件存在
list.files("01_raw_data/02_outcome")
# ============================================================
# 4. FinnGen outcome VCF
# ============================================================

file_finngen <-
  "01_raw_data/02_outcome/finn-b-C3_STOMACH.vcf.gz"

file.exists(file_finngen)
# ============================================================
# 5. Extract prespecified IVs from FinnGen
# ============================================================

finngen_result <- extract_outcome_snps(
  vcf_file = file_finngen,
  target_snps = target_snps_finn,
  chunk_size = 50000
)

finngen_raw <- finngen_result$data


# Check matching
finngen_result$n_requested
finngen_result$n_found
finngen_result$n_missing
finngen_result$missing
# ============================================================
# 7. FinnGen outcome data QC
# ============================================================

nrow(finngen_raw)

# Check missing values
sum(is.na(finngen_raw$beta))
sum(is.na(finngen_raw$se))
sum(is.na(finngen_raw$pval))
sum(is.na(finngen_raw$eaf))

# Check duplicated SNPs
sum(duplicated(finngen_raw$SNP))

# Inspect extracted outcome data
finngen_raw
# ============================================================
# 8. Save extracted FinnGen outcome data
# ============================================================

saveRDS(
  finngen_raw,
  "02_intermediate_data/02_outcome/FinnGen/finn-b-C3_STOMACH_extracted_IVs.rds"
)

write.csv(
  finngen_raw,
  "02_intermediate_data/02_outcome/FinnGen/finn-b-C3_STOMACH_extracted_IVs.csv",
  row.names = FALSE
)
# ============================================================
# 9. Format FinnGen outcome for TwoSampleMR
# ============================================================

finngen_format <- finngen_raw

finngen_format$Phenotype <- "FinnGen gastric cancer"

finngen_outcome <- TwoSampleMR::format_data(
  finngen_format,
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

nrow(finngen_outcome)
names(finngen_outcome)
# ============================================================
# 10. Separate FinnGen outcome by exposure
# ============================================================

finngen_veillonella <- finngen_outcome[
  finngen_outcome$SNP %in% veillonella_iv$SNP,
]

finngen_veillonellaceae <- finngen_outcome[
  finngen_outcome$SNP %in% veillonellaceae_iv$SNP,
]

nrow(finngen_veillonella)
nrow(finngen_veillonellaceae)
#preregistration 已经冻结为保守处理 ambiguous palindromic SNP，仍然使用
action = 3
# ============================================================
# 11. Harmonisation
# ============================================================

harm_veillonella_finn <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonella_iv,
  outcome_dat = finngen_veillonella,
  action = 3
)

harm_veillonellaceae_finn <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonellaceae_iv,
  outcome_dat = finngen_veillonellaceae,
  action = 3
)
# ============================================================
# 12. Harmonisation QC
# ============================================================

table(harm_veillonella_finn$mr_keep)
table(harm_veillonellaceae_finn$mr_keep)

harm_veillonella_finn[
  harm_veillonella_finn$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]

harm_veillonellaceae_finn[
  harm_veillonellaceae_finn$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]
#生成MR-ready数据
veillonella_finn_ready <- subset(
  harm_veillonella_finn,
  mr_keep == TRUE
)

veillonellaceae_finn_ready <- subset(
  harm_veillonellaceae_finn,
  mr_keep == TRUE
)

nrow(veillonella_finn_ready)
nrow(veillonellaceae_finn_ready)

sum(is.na(veillonella_finn_ready$beta.exposure))
sum(is.na(veillonella_finn_ready$beta.outcome))

sum(is.na(veillonellaceae_finn_ready$beta.exposure))
sum(is.na(veillonellaceae_finn_ready$beta.outcome))
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
# 14. FinnGen MR: Veillonella -> gastric cancer
# ============================================================

mr_veillonella_finn <- TwoSampleMR::mr(
  veillonella_finn_ready,
  method_list = mr_methods
)

mr_veillonella_finn


# ============================================================
# 15. FinnGen MR: Veillonellaceae -> gastric cancer
# ============================================================

mr_veillonellaceae_finn <- TwoSampleMR::mr(
  veillonellaceae_finn_ready,
  method_list = mr_methods
)

mr_veillonellaceae_finn
#转换OR值
# ============================================================
# 16. Convert to odds ratios
# ============================================================

mr_veillonella_finn_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonella_finn
  )

mr_veillonellaceae_finn_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonellaceae_finn
  )

mr_veillonella_finn_or
mr_veillonellaceae_finn_or
# ============================================================
# 17. Heterogeneity
# ============================================================

het_veillonella_finn <- TwoSampleMR::mr_heterogeneity(
  veillonella_finn_ready
)

het_veillonellaceae_finn <- TwoSampleMR::mr_heterogeneity(
  veillonellaceae_finn_ready
)

het_veillonella_finn
het_veillonellaceae_finn


# ============================================================
# 18. MR-Egger intercept
# ============================================================

pleio_veillonella_finn <- TwoSampleMR::mr_pleiotropy_test(
  veillonella_finn_ready
)

pleio_veillonellaceae_finn <- TwoSampleMR::mr_pleiotropy_test(
  veillonellaceae_finn_ready
)

pleio_veillonella_finn
pleio_veillonellaceae_finn
# ============================================================
# 19. Leave-one-out
# ============================================================

loo_veillonella_finn <- TwoSampleMR::mr_leaveoneout(
  veillonella_finn_ready
)

loo_veillonellaceae_finn <- TwoSampleMR::mr_leaveoneout(
  veillonellaceae_finn_ready
)

loo_veillonella_finn
loo_veillonellaceae_finn


# ============================================================
# 20. Single-SNP analysis
# ============================================================

single_veillonella_finn <- TwoSampleMR::mr_singlesnp(
  veillonella_finn_ready
)

single_veillonellaceae_finn <- TwoSampleMR::mr_singlesnp(
  veillonellaceae_finn_ready
)

single_veillonella_finn
single_veillonellaceae_finn
# ============================================================
# 21. MR-PRESSO
# ============================================================

library(MRPRESSO)

presso_veillonella_finn <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonella_finn_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonellaceae_finn <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonellaceae_finn_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonella_finn
presso_veillonellaceae_finn
# ============================================================
# 22. Save FinnGen MR results
# ============================================================

write.csv(
  mr_veillonella_finn_or,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_MR.csv",
  row.names = FALSE
)

write.csv(
  mr_veillonellaceae_finn_or,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_MR.csv",
  row.names = FALSE
)

write.csv(
  het_veillonella_finn,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  het_veillonellaceae_finn,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonella_finn,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonellaceae_finn,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonella_finn,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonellaceae_finn,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  single_veillonella_finn,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  single_veillonellaceae_finn,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_singleSNP.csv",
  row.names = FALSE
)

saveRDS(
  presso_veillonella_finn,
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_MRPRESSO.rds"
)

saveRDS(
  presso_veillonellaceae_finn,
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_MRPRESSO.rds"
)

capture.output(
  presso_veillonella_finn,
  file = "05_results/03_finngen_sensitivity/Veillonella_FinnGen_MRPRESSO.txt"
)

capture.output(
  presso_veillonellaceae_finn,
  file = "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_MRPRESSO.txt"
)
# ============================================================
# 23. Save FinnGen figures
# ============================================================

forest_veillonella_finn <-
  TwoSampleMR::mr_forest_plot(single_veillonella_finn)

forest_veillonellaceae_finn <-
  TwoSampleMR::mr_forest_plot(single_veillonellaceae_finn)

loo_plot_veillonella_finn <-
  TwoSampleMR::mr_leaveoneout_plot(loo_veillonella_finn)

loo_plot_veillonellaceae_finn <-
  TwoSampleMR::mr_leaveoneout_plot(loo_veillonellaceae_finn)

library(ggplot2)
ggsave(
  "06_figures/03_finngen_sensitivity/Veillonella_FinnGen_singleSNP_forest.pdf",
  forest_veillonella_finn[[1]],
  width = 7,
  height = 6
)

ggsave(
  "06_figures/03_finngen_sensitivity/Veillonellaceae_FinnGen_singleSNP_forest.pdf",
  forest_veillonellaceae_finn[[1]],
  width = 7,
  height = 8
)

ggsave(
  "06_figures/03_finngen_sensitivity/Veillonella_FinnGen_leaveoneout.pdf",
  loo_plot_veillonella_finn[[1]],
  width = 7,
  height = 6
)

ggsave(
  "06_figures/03_finngen_sensitivity/Veillonellaceae_FinnGen_leaveoneout.pdf",
  loo_plot_veillonellaceae_finn[[1]],
  width = 7,
  height = 8
)
#检查保存项目
list.files("06_figures/03_finngen_sensitivity")
list.files("05_results/03_finngen_sensitivity")
