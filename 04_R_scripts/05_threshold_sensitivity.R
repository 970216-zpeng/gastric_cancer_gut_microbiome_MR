# ============================================================
# 05_threshold_sensitivity.R
# Project: 01_MR_GastricCancer
# Sensitivity threshold: P <= 1e-5
# Primary outcome: GCST90018849
# ============================================================

source("04_R_scripts/00_functions.R")

library(TwoSampleMR)


# ============================================================
# 1. Load sensitivity exposure IVs
# ============================================================

veillonella_sens_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonella_IV_sensitivity_p1e5.rds"
)

veillonellaceae_sens_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_sensitivity_p1e5.rds"
)


nrow(veillonella_sens_iv)
nrow(veillonellaceae_sens_iv)
# ============================================================
# 2. Create target SNP list
# ============================================================

target_snps_sens <- unique(
  c(
    veillonella_sens_iv$SNP,
    veillonellaceae_sens_iv$SNP
  )
)

length(target_snps_sens)

target_snps_sens
#检查
intersect(
  veillonella_sens_iv$SNP,
  veillonellaceae_sens_iv$SNP
)

intersect(
  veillonella_sens_iv$SNP,
  veillonellaceae_sens_iv$SNP
)
# ============================================================
# 3. Extract sensitivity IVs from primary GC outcome
# ============================================================

file_gc_primary <-
  "01_raw_data/02_outcome/ebi-a-GCST90018849.vcf.gz"

gc_sens_result <- extract_outcome_snps(
  vcf_file = file_gc_primary,
  target_snps = target_snps_sens,
  chunk_size = 50000
)

gc_sens <- gc_sens_result$data
#检查匹配情况
gc_sens_result$n_requested
gc_sens_result$n_found
gc_sens_result$n_missing
gc_sens_result$missing
#在检查完整性
sum(is.na(gc_sens$beta))
sum(is.na(gc_sens$se))
sum(is.na(gc_sens$pval))
sum(is.na(gc_sens$eaf))
sum(duplicated(gc_sens$SNP))
# ============================================================
# 4. Format sensitivity outcome
# ============================================================

gc_sens$Phenotype <- "Gastric cancer"

gc_sens_outcome <- TwoSampleMR::format_data(
  gc_sens,
  type = "outcome",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "pos"
)

gc_sens_outcome$id.outcome <- "GCST90018849"
gc_sens_outcome$outcome <- "Gastric cancer"
# ============================================================
# 5. Identify palindromic SNPs
# ============================================================

is_palindromic <- function(a1, a2) {
  pair <- paste0(
    pmin(a1, a2),
    pmax(a1, a2)
  )
  pair %in% c("AT", "CG")
}

veillonella_sens_iv$palindromic_check <- is_palindromic(
  veillonella_sens_iv$effect_allele.exposure,
  veillonella_sens_iv$other_allele.exposure
)

veillonellaceae_sens_iv$palindromic_check <- is_palindromic(
  veillonellaceae_sens_iv$effect_allele.exposure,
  veillonellaceae_sens_iv$other_allele.exposure
)

sum(veillonella_sens_iv$palindromic_check)
sum(veillonellaceae_sens_iv$palindromic_check)
#具体SNP
veillonella_sens_iv[
  veillonella_sens_iv$palindromic_check,
  c("SNP",
    "effect_allele.exposure",
    "other_allele.exposure")
]

veillonellaceae_sens_iv[
  veillonellaceae_sens_iv$palindromic_check,
  c("SNP",
    "effect_allele.exposure",
    "other_allele.exposure")
]
# ============================================================
# 6. Harmonisation
# ============================================================

harm_veillonella_sens <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonella_sens_iv,
  outcome_dat = gc_sens_outcome,
  action = 3
)

harm_veillonellaceae_sens <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonellaceae_sens_iv,
  outcome_dat = gc_sens_outcome,
  action = 3
)
#看保留多少
table(harm_veillonella_sens$mr_keep)
table(harm_veillonellaceae_sens$mr_keep)
#查看删除的SNP
harm_veillonella_sens[
  harm_veillonella_sens$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]

harm_veillonellaceae_sens[
  harm_veillonellaceae_sens$mr_keep == FALSE,
  c("SNP", "palindromic", "ambiguous", "mr_keep")
]
# ============================================================
# 7. Create sensitivity MR-ready datasets
# ============================================================

veillonella_sens_ready <- subset(
  harm_veillonella_sens,
  mr_keep == TRUE
)

veillonellaceae_sens_ready <- subset(
  harm_veillonellaceae_sens,
  mr_keep == TRUE
)

nrow(veillonella_sens_ready)
nrow(veillonellaceae_sens_ready)
#完整性检查
sum(is.na(veillonella_sens_ready$beta.exposure))
sum(is.na(veillonella_sens_ready$beta.outcome))

sum(is.na(veillonellaceae_sens_ready$beta.exposure))
sum(is.na(veillonellaceae_sens_ready$beta.outcome))
#保存完整harmonisation audit和MR-ready数据
dir.create(
  "02_intermediate_data/04_threshold_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  harm_veillonella_sens,
  "02_intermediate_data/04_threshold_sensitivity/Veillonella_p1e5_GCST90018849_harmonisation_audit.csv",
  row.names = FALSE
)

write.csv(
  harm_veillonellaceae_sens,
  "02_intermediate_data/04_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_harmonisation_audit.csv",
  row.names = FALSE
)

saveRDS(
  veillonella_sens_ready,
  "02_intermediate_data/04_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MR_ready.rds"
)

saveRDS(
  veillonellaceae_sens_ready,
  "02_intermediate_data/04_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MR_ready.rds"
)
# ============================================================
# 8. Threshold sensitivity MR
# ============================================================

mr_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)

mr_veillonella_sens <- TwoSampleMR::mr(
  veillonella_sens_ready,
  method_list = mr_methods
)

mr_veillonellaceae_sens <- TwoSampleMR::mr(
  veillonellaceae_sens_ready,
  method_list = mr_methods
)

mr_veillonella_sens
mr_veillonellaceae_sens
#转换or
mr_veillonella_sens_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonella_sens
  )

mr_veillonellaceae_sens_or <-
  TwoSampleMR::generate_odds_ratios(
    mr_veillonellaceae_sens
  )

mr_veillonella_sens_or
mr_veillonellaceae_sens_or
# ============================================================
# 9. Threshold sensitivity heterogeneity
# ============================================================

het_veillonella_sens <- TwoSampleMR::mr_heterogeneity(
  veillonella_sens_ready
)

het_veillonellaceae_sens <- TwoSampleMR::mr_heterogeneity(
  veillonellaceae_sens_ready
)

het_veillonella_sens
het_veillonellaceae_sens
# ============================================================
# 10. Threshold sensitivity pleiotropy
# ============================================================

pleio_veillonella_sens <- TwoSampleMR::mr_pleiotropy_test(
  veillonella_sens_ready
)

pleio_veillonellaceae_sens <- TwoSampleMR::mr_pleiotropy_test(
  veillonellaceae_sens_ready
)

pleio_veillonella_sens
pleio_veillonellaceae_sens
# ============================================================
# 11. Leave-one-out analysis
# ============================================================

loo_veillonella_sens <- TwoSampleMR::mr_leaveoneout(
  veillonella_sens_ready
)

loo_veillonellaceae_sens <- TwoSampleMR::mr_leaveoneout(
  veillonellaceae_sens_ready
)

loo_veillonella_sens
loo_veillonellaceae_sens
# ============================================================
# 12. Single-SNP analysis
# ============================================================

single_veillonella_sens <- TwoSampleMR::mr_singlesnp(
  veillonella_sens_ready
)

single_veillonellaceae_sens <- TwoSampleMR::mr_singlesnp(
  veillonellaceae_sens_ready
)

single_veillonella_sens
single_veillonellaceae_sens
# ============================================================
# 13. Sensitivity QC plots
# ============================================================

forest_veillonella_sens <-
  TwoSampleMR::mr_forest_plot(
    single_veillonella_sens
  )

forest_veillonellaceae_sens <-
  TwoSampleMR::mr_forest_plot(
    single_veillonellaceae_sens
  )


loo_plot_veillonella_sens <-
  TwoSampleMR::mr_leaveoneout_plot(
    loo_veillonella_sens
  )

loo_plot_veillonellaceae_sens <-
  TwoSampleMR::mr_leaveoneout_plot(
    loo_veillonellaceae_sens
  )
#建立 sensitivity results / figures 文件夹
dir.create(
  "05_results/02_threshold_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "06_figures/02_threshold_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)
#保存数值结果
write.csv(
  mr_veillonella_sens_or,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MR.csv",
  row.names = FALSE
)

write.csv(
  mr_veillonellaceae_sens_or,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MR.csv",
  row.names = FALSE
)

write.csv(
  het_veillonella_sens,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  het_veillonellaceae_sens,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonella_sens,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonellaceae_sens,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonella_sens,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonellaceae_sens,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  single_veillonella_sens,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  single_veillonellaceae_sens,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_singleSNP.csv",
  row.names = FALSE
)
#保存图片forest / leave-one-out
ggplot2::ggsave(
  "06_figures/02_threshold_sensitivity/Veillonella_p1e5_singleSNP_forest.pdf",
  plot = forest_veillonella_sens[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/02_threshold_sensitivity/Veillonellaceae_p1e5_singleSNP_forest.pdf",
  plot = forest_veillonellaceae_sens[[1]],
  width = 7,
  height = 8
)

ggplot2::ggsave(
  "06_figures/02_threshold_sensitivity/Veillonella_p1e5_leaveoneout.pdf",
  plot = loo_plot_veillonella_sens[[1]],
  width = 7,
  height = 6
)

ggplot2::ggsave(
  "06_figures/02_threshold_sensitivity/Veillonellaceae_p1e5_leaveoneout.pdf",
  plot = loo_plot_veillonellaceae_sens[[1]],
  width = 7,
  height = 8
)
# ============================================================
# 14. MR-PRESSO
# ============================================================

library(MRPRESSO)

presso_veillonellaceae_sens <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonellaceae_sens_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonellaceae_sens
#先保存
saveRDS(
  presso_veillonellaceae_sens,
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MRPRESSO.rds"
)

capture.output(
  presso_veillonellaceae_sens,
  file =
    "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MRPRESSO.txt"
)
# ============================================================
# 15. MR-PRESSO: Veillonella threshold sensitivity
# ============================================================

presso_veillonella_sens <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonella_sens_ready,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonella_sens
#保存
saveRDS(
  presso_veillonella_sens,
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MRPRESSO.rds"
)

capture.output(
  presso_veillonella_sens,
  file =
    "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MRPRESSO.txt"
)
#检查
list.files("05_results/02_threshold_sensitivity")
list.files("06_figures/02_threshold_sensitivity")
list.files("05_results/01_primary_MR")
list.files("06_figures/01_primary_MR")


