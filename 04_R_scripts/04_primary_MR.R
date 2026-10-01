# ============================================================
# 04_primary_MR.R
# Project: 01_MR_GastricCancer
# Primary outcome: GCST90018849
# Purpose:
#   Primary MR analysis
# ============================================================

library(TwoSampleMR)


# ============================================================
# 1. Load MR-ready datasets
# ============================================================

veillonella_dat <- readRDS(
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_MR_ready.rds"
)

veillonellaceae_dat <- readRDS(
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_MR_ready.rds"
)


# Check
nrow(veillonella_dat)
nrow(veillonellaceae_dat)
# ============================================================
# 2. Define MR methods
# ============================================================

mr_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)
# ============================================================
# 3. MR: Veillonella -> Gastric cancer
# ============================================================

mr_veillonella <- TwoSampleMR::mr(
  veillonella_dat,
  method_list = mr_methods
)

mr_veillonella
# ============================================================
# 4. MR: Veillonellaceae -> Gastric cancer
# ============================================================

mr_veillonellaceae <- TwoSampleMR::mr(
  veillonellaceae_dat,
  method_list = mr_methods
)

mr_veillonellaceae
# ============================================================
# 5. Convert MR estimates to odds ratios
# ============================================================

mr_veillonella_or <- TwoSampleMR::generate_odds_ratios(
  mr_veillonella
)

mr_veillonellaceae_or <- TwoSampleMR::generate_odds_ratios(
  mr_veillonellaceae
)

mr_veillonella_or
mr_veillonellaceae_or
# ============================================================
# 6. Heterogeneity
# ============================================================

het_veillonella <- TwoSampleMR::mr_heterogeneity(
  veillonella_dat
)

het_veillonellaceae <- TwoSampleMR::mr_heterogeneity(
  veillonellaceae_dat
)

het_veillonella
het_veillonellaceae
# ============================================================
# 7. Horizontal pleiotropy
# ============================================================

pleio_veillonella <- TwoSampleMR::mr_pleiotropy_test(
  veillonella_dat
)

pleio_veillonellaceae <- TwoSampleMR::mr_pleiotropy_test(
  veillonellaceae_dat
)

pleio_veillonella
pleio_veillonellaceae
#检查
het_veillonella
het_veillonellaceae
pleio_veillonella
pleio_veillonellaceae
# ============================================================
# 9. Leave-one-out analysis
# ============================================================

loo_veillonella <- TwoSampleMR::mr_leaveoneout(
  veillonella_dat
)

loo_veillonellaceae <- TwoSampleMR::mr_leaveoneout(
  veillonellaceae_dat
)

loo_veillonella
loo_veillonellaceae
# ============================================================
# 10. Single-SNP analysis
# ============================================================

single_veillonella <- TwoSampleMR::mr_singlesnp(
  veillonella_dat
)

single_veillonellaceae <- TwoSampleMR::mr_singlesnp(
  veillonellaceae_dat
)

single_veillonella

single_veillonellaceae
# ============================================================
# 11. Single-SNP forest plots
# ============================================================

forest_veillonella <-
  TwoSampleMR::mr_forest_plot(
    single_veillonella
  )

forest_veillonellaceae <-
  TwoSampleMR::mr_forest_plot(
    single_veillonellaceae
  )

forest_veillonella[[1]]

forest_veillonellaceae[[1]]
# ============================================================
# 12. Leave-one-out plots
# ============================================================

loo_plot_veillonella <-
  TwoSampleMR::mr_leaveoneout_plot(
    loo_veillonella
  )

loo_plot_veillonellaceae <-
  TwoSampleMR::mr_leaveoneout_plot(
    loo_veillonellaceae
  )

loo_plot_veillonella[[1]]

loo_plot_veillonellaceae[[1]]
# ============================================================
# Check plotting warnings
# ============================================================

single_veillonella[
  !complete.cases(single_veillonella[, c("b", "se")]),
]

single_veillonellaceae[
  !complete.cases(single_veillonellaceae[, c("b", "se")]),
]

loo_veillonella[
  !complete.cases(loo_veillonella[, c("b", "se")]),
]

loo_veillonellaceae[
  !complete.cases(loo_veillonellaceae[, c("b", "se")]),
]
#再检查所有估计值是不是 finite
all(is.finite(single_veillonella$b))
all(is.finite(single_veillonella$se))

all(is.finite(single_veillonellaceae$b))
all(is.finite(single_veillonellaceae$se))

all(is.finite(loo_veillonella$b))
all(is.finite(loo_veillonella$se))

all(is.finite(loo_veillonellaceae$b))
all(is.finite(loo_veillonellaceae$se))
#确认
nrow(single_veillonella)
nrow(single_veillonellaceae)

nrow(loo_veillonella)
nrow(loo_veillonellaceae)
# ============================================================
# 13. Save primary MR sensitivity figures
# ============================================================

dir.create(
  "06_figures/01_primary_MR",
  recursive = TRUE,
  showWarnings = FALSE
)
#保存 Single-SNP forest plots
ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonella_GCST90018849_singleSNP_forest.pdf",
  plot = forest_veillonella[[1]],
  width = 7,
  height = 5
)

ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonellaceae_GCST90018849_singleSNP_forest.pdf",
  plot = forest_veillonellaceae[[1]],
  width = 7,
  height = 6
)
#同时保存 300 dpi PNG
ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonella_GCST90018849_singleSNP_forest.png",
  plot = forest_veillonella[[1]],
  width = 7,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonellaceae_GCST90018849_singleSNP_forest.png",
  plot = forest_veillonellaceae[[1]],
  width = 7,
  height = 6,
  dpi = 300
)
#保存 Leave-one-out plots
ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonella_GCST90018849_leaveoneout.pdf",
  plot = loo_plot_veillonella[[1]],
  width = 7,
  height = 5
)

ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonellaceae_GCST90018849_leaveoneout.pdf",
  plot = loo_plot_veillonellaceae[[1]],
  width = 7,
  height = 6
)
#PNG
ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonella_GCST90018849_leaveoneout.png",
  plot = loo_plot_veillonella[[1]],
  width = 7,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  "06_figures/01_primary_MR/Veillonellaceae_GCST90018849_leaveoneout.png",
  plot = loo_plot_veillonellaceae[[1]],
  width = 7,
  height = 6,
  dpi = 300
)

#对应的原始数据也保存
dir.create(
  "05_results/01_primary_MR",
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(
  single_veillonella,
  "05_results/01_primary_MR/Veillonella_GCST90018849_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  single_veillonellaceae,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_singleSNP.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonella,
  "05_results/01_primary_MR/Veillonella_GCST90018849_leaveoneout.csv",
  row.names = FALSE
)

write.csv(
  loo_veillonellaceae,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_leaveoneout.csv",
  row.names = FALSE
)
#确认
list.files("05_results/01_primary_MR")

# ============================================================
# 14. MR-PRESSO
# ============================================================

# Check package
requireNamespace("MRPRESSO", quietly = TRUE)
print(requireNamespace("MRPRESSO", quietly = TRUE))
library(MRPRESSO)
#
presso_veillonellaceae <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonellaceae_dat,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonellaceae
#保存MR-PRESSO
saveRDS(
  presso_veillonellaceae,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_MRPRESSO.rds"
)

capture.output(
  presso_veillonellaceae,
  file =
    "05_results/01_primary_MR/Veillonellaceae_GCST90018849_MRPRESSO.txt"
)
#
presso_veillonella <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = veillonella_dat,
  NbDistribution = 10000,
  SignifThreshold = 0.05
)

presso_veillonella
#保存MR-PRESSO
saveRDS(
  presso_veillonella,
  "05_results/01_primary_MR/Veillonella_GCST90018849_MRPRESSO.rds"
)

capture.output(
  presso_veillonella,
  file =
    "05_results/01_primary_MR/Veillonella_GCST90018849_MRPRESSO.txt"
)
# ============================================================
# Save primary MR summary results
# ============================================================

write.csv(
  mr_veillonella_or,
  "05_results/01_primary_MR/Veillonella_GCST90018849_MR.csv",
  row.names = FALSE
)

write.csv(
  mr_veillonellaceae_or,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_MR.csv",
  row.names = FALSE
)

write.csv(
  het_veillonella,
  "05_results/01_primary_MR/Veillonella_GCST90018849_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  het_veillonellaceae,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_heterogeneity.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonella,
  "05_results/01_primary_MR/Veillonella_GCST90018849_Egger_intercept.csv",
  row.names = FALSE
)

write.csv(
  pleio_veillonellaceae,
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_Egger_intercept.csv",
  row.names = FALSE
)
#检查
list.files("05_results/01_primary_MR")
