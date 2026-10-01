# ============================================================
# 09_summary_results.R
# Project: 01_MR_GastricCancer
#
# Purpose:
#   Summarise all preregistered MR analyses
#
# Includes:
#   1. Primary forward MR
#   2. Alternative-threshold sensitivity
#   3. FinnGen sensitivity
#   4. East Asian cross-ancestry sensitivity
#   5. Reverse MR
# ============================================================

library(dplyr)
library(ggplot2)

dir.create(
  "05_results/06_summary",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "06_figures/06_summary",
  recursive = TRUE,
  showWarnings = FALSE
)
# ============================================================
# 1. Helper function: read MR results
# ============================================================

read_mr_result <- function(
    file,
    analysis,
    exposure_name,
    outcome_name,
    ancestry,
    direction,
    analysis_role
) {
  
  x <- read.csv(file, stringsAsFactors = FALSE)
  
  x$analysis <- analysis
  x$exposure_label <- exposure_name
  x$outcome_label <- outcome_name
  x$ancestry <- ancestry
  x$direction <- direction
  x$analysis_role <- analysis_role
  
  x
}
# ============================================================
# 2. Primary forward MR
# ============================================================

primary_veillonella <- read_mr_result(
  "05_results/01_primary_MR/Veillonella_GCST90018849_MR.csv",
  analysis = "Primary GCST90018849",
  exposure_name = "Veillonella",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Primary"
)

primary_veillonellaceae <- read_mr_result(
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_MR.csv",
  analysis = "Primary GCST90018849",
  exposure_name = "Veillonellaceae",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Primary"
)
# ============================================================
# 3. Alternative instrument-threshold sensitivity
# ============================================================

threshold_veillonella <- read_mr_result(
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MR.csv",
  analysis = "Threshold sensitivity P<=1e-5",
  exposure_name = "Veillonella",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Sensitivity"
)

threshold_veillonellaceae <- read_mr_result(
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MR.csv",
  analysis = "Threshold sensitivity P<=1e-5",
  exposure_name = "Veillonellaceae",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Sensitivity"
)
# ============================================================
# 4. Independent European FinnGen sensitivity
# ============================================================

finngen_veillonella <- read_mr_result(
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_MR.csv",
  analysis = "FinnGen",
  exposure_name = "Veillonella",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Sensitivity"
)

finngen_veillonellaceae <- read_mr_result(
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_MR.csv",
  analysis = "FinnGen",
  exposure_name = "Veillonellaceae",
  outcome_name = "Gastric cancer",
  ancestry = "European",
  direction = "Forward",
  analysis_role = "Sensitivity"
)
# ============================================================
# 5. East Asian cross-ancestry sensitivity
# ============================================================

eas_veillonella <- read_mr_result(
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_MR.csv",
  analysis = "East Asian GCST90018629",
  exposure_name = "Veillonella",
  outcome_name = "Gastric cancer",
  ancestry = "East Asian",
  direction = "Forward",
  analysis_role = "Sensitivity"
)

eas_veillonellaceae <- read_mr_result(
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_MR.csv",
  analysis = "East Asian GCST90018629",
  exposure_name = "Veillonellaceae",
  outcome_name = "Gastric cancer",
  ancestry = "East Asian",
  direction = "Forward",
  analysis_role = "Sensitivity"
)
# ============================================================
# 6. Reverse MR
# ============================================================

reverse_veillonella <- read_mr_result(
  "05_results/05_reverse_MR/GC_to_Veillonella_MR.csv",
  analysis = "Reverse MR",
  exposure_name = "Gastric cancer",
  outcome_name = "Veillonella",
  ancestry = "European",
  direction = "Reverse",
  analysis_role = "Reverse"
)

reverse_veillonellaceae <- read_mr_result(
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_MR.csv",
  analysis = "Reverse MR",
  exposure_name = "Gastric cancer",
  outcome_name = "Veillonellaceae",
  ancestry = "European",
  direction = "Reverse",
  analysis_role = "Reverse"
)
# ============================================================
# 7. Combine all MR results
# ============================================================

all_mr <- bind_rows(
  primary_veillonella,
  primary_veillonellaceae,
  threshold_veillonella,
  threshold_veillonellaceae,
  finngen_veillonella,
  finngen_veillonellaceae,
  eas_veillonella,
  eas_veillonellaceae,
  reverse_veillonella,
  reverse_veillonellaceae
)

dim(all_mr)
table(all_mr$analysis)
table(all_mr$direction)
#保存
write.csv(
  all_mr,
  "05_results/06_summary/all_MR_results.csv",
  row.names = FALSE
)
# ============================================================
# 8. IVW summary
# ============================================================

ivw_summary <- all_mr %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    significance_threshold = 0.025,
    significant = pval < significance_threshold
  )

ivw_summary
#对于 reverse MR没有 OR，所以我们补一个统一的展示字段
ivw_summary <- ivw_summary %>%
  mutate(
    estimate_type = ifelse(
      direction == "Forward",
      "OR",
      "Beta"
    ),
    
    estimate = ifelse(
      direction == "Forward",
      or,
      b
    ),
    
    ci_lower = ifelse(
      direction == "Forward",
      or_lci95,
      b - 1.96 * se
    ),
    
    ci_upper = ifelse(
      direction == "Forward",
      or_uci95,
      b + 1.96 * se
    )
  )
#保存需要列
ivw_table <- ivw_summary %>%
  select(
    analysis,
    analysis_role,
    direction,
    ancestry,
    exposure_label,
    outcome_label,
    nsnp,
    estimate_type,
    estimate,
    ci_lower,
    ci_upper,
    pval,
    significant
  )

ivw_table
#保存
write.csv(
  ivw_table,
  "05_results/06_summary/IVW_summary.csv",
  row.names = FALSE
)
#增加一个“论文可直接阅读”的格式化列
ivw_table_formatted <- ivw_table %>%
  mutate(
    estimate_95CI = sprintf(
      "%.3f (%.3f–%.3f)",
      estimate,
      ci_lower,
      ci_upper
    ),

    p_display = ifelse(
      pval < 0.001,
      "<0.001",
      sprintf("%.3f", pval)
    )
  )

ivw_table_formatted
#保存
write.csv(
  ivw_table_formatted,
  "05_results/06_summary/IVW_summary_formatted.csv",
  row.names = FALSE
)

# ============================================================
# 9. Direction consistency
# ============================================================

veillonellaceae_forward <- ivw_summary %>%
  filter(
    direction == "Forward",
    exposure_label == "Veillonellaceae"
  ) %>%
  select(
    analysis,
    b,
    or,
    or_lci95,
    or_uci95,
    pval
  )

veillonellaceae_forward
#在检查all
all(veillonellaceae_forward$b > 0)
# ============================================================
# 10. Helper function: summarise QC results
# ============================================================

read_qc <- function(
    heterogeneity_file,
    egger_file,
    presso_file,
    analysis,
    comparison
) {
  
  het <- read.csv(
    heterogeneity_file,
    stringsAsFactors = FALSE
  )
  
  egger <- read.csv(
    egger_file,
    stringsAsFactors = FALSE
  )
  
  presso <- readRDS(presso_file)
  
  # IVW heterogeneity
  ivw_q <- het[
    het$method == "Inverse variance weighted",
  ]
  
  # MR-PRESSO global-test P value
  presso_global_p <-
    presso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  
  data.frame(
    analysis = analysis,
    comparison = comparison,
    
    IVW_Q = ivw_q$Q[1],
    IVW_Q_df = ivw_q$Q_df[1],
    IVW_Q_p = ivw_q$Q_pval[1],
    
    Egger_intercept = egger$egger_intercept[1],
    Egger_intercept_SE = egger$se[1],
    Egger_intercept_p = egger$pval[1],
    
    MRPRESSO_global_p = as.numeric(presso_global_p),
    
    stringsAsFactors = FALSE
  )
}
#测试上面函数
qc_primary_veillonella <- read_qc(
  "05_results/01_primary_MR/Veillonella_GCST90018849_heterogeneity.csv",
  "05_results/01_primary_MR/Veillonella_GCST90018849_Egger_intercept.csv",
  "05_results/01_primary_MR/Veillonella_GCST90018849_MRPRESSO.rds",
  analysis = "Primary GCST90018849",
  comparison = "Veillonella -> Gastric cancer"
)

qc_primary_veillonellaceae <- read_qc(
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_heterogeneity.csv",
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_Egger_intercept.csv",
  "05_results/01_primary_MR/Veillonellaceae_GCST90018849_MRPRESSO.rds",
  analysis = "Primary GCST90018849",
  comparison = "Veillonellaceae -> Gastric cancer"
)

qc_primary_veillonella
qc_primary_veillonellaceae
# ============================================================
# 11. Threshold sensitivity QC
# ============================================================

qc_threshold_veillonella <- read_qc(
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_heterogeneity.csv",
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_Egger_intercept.csv",
  "05_results/02_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MRPRESSO.rds",
  analysis = "Threshold sensitivity P<=1e-5",
  comparison = "Veillonella -> Gastric cancer"
)

qc_threshold_veillonellaceae <- read_qc(
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_heterogeneity.csv",
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_Egger_intercept.csv",
  "05_results/02_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MRPRESSO.rds",
  analysis = "Threshold sensitivity P<=1e-5",
  comparison = "Veillonellaceae -> Gastric cancer"
)
# ============================================================
# 12. FinnGen QC
# ============================================================

qc_finngen_veillonella <- read_qc(
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_heterogeneity.csv",
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_Egger_intercept.csv",
  "05_results/03_finngen_sensitivity/Veillonella_FinnGen_MRPRESSO.rds",
  analysis = "FinnGen",
  comparison = "Veillonella -> Gastric cancer"
)

qc_finngen_veillonellaceae <- read_qc(
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_heterogeneity.csv",
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_Egger_intercept.csv",
  "05_results/03_finngen_sensitivity/Veillonellaceae_FinnGen_MRPRESSO.rds",
  analysis = "FinnGen",
  comparison = "Veillonellaceae -> Gastric cancer"
)
# ============================================================
# 13. East Asian QC
# ============================================================

qc_eas_veillonella <- read_qc(
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_heterogeneity.csv",
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_Egger_intercept.csv",
  "05_results/04_east_asian_sensitivity/Veillonella_GCST90018629_MRPRESSO.rds",
  analysis = "East Asian GCST90018629",
  comparison = "Veillonella -> Gastric cancer"
)

qc_eas_veillonellaceae <- read_qc(
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_heterogeneity.csv",
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_Egger_intercept.csv",
  "05_results/04_east_asian_sensitivity/Veillonellaceae_GCST90018629_MRPRESSO.rds",
  analysis = "East Asian GCST90018629",
  comparison = "Veillonellaceae -> Gastric cancer"
)
# ============================================================
# 14. Reverse MR QC
# ============================================================

qc_reverse_veillonella <- read_qc(
  "05_results/05_reverse_MR/GC_to_Veillonella_heterogeneity.csv",
  "05_results/05_reverse_MR/GC_to_Veillonella_Egger_intercept.csv",
  "05_results/05_reverse_MR/GC_to_Veillonella_MRPRESSO.rds",
  analysis = "Reverse MR",
  comparison = "Gastric cancer -> Veillonella"
)

qc_reverse_veillonellaceae <- read_qc(
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_heterogeneity.csv",
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_Egger_intercept.csv",
  "05_results/05_reverse_MR/GC_to_Veillonellaceae_MRPRESSO.rds",
  analysis = "Reverse MR",
  comparison = "Gastric cancer -> Veillonellaceae"
)
# ============================================================
# 15. Combine QC results
# ============================================================

qc_summary <- dplyr::bind_rows(
  qc_primary_veillonella,
  qc_primary_veillonellaceae,
  qc_threshold_veillonella,
  qc_threshold_veillonellaceae,
  qc_finngen_veillonella,
  qc_finngen_veillonellaceae,
  qc_eas_veillonella,
  qc_eas_veillonellaceae,
  qc_reverse_veillonella,
  qc_reverse_veillonellaceae
)

qc_summary
#检查
nrow(qc_summary)
#增加 3 个自动判断字段：这样以后写论文不需要人工逐个查 P 值
qc_summary <- qc_summary %>%
  mutate(
    heterogeneity_evidence = IVW_Q_p < 0.05,
    directional_pleiotropy_evidence = Egger_intercept_p < 0.05,
    presso_global_evidence = MRPRESSO_global_p < 0.05
  )

qc_summary
#检查
table(qc_summary$heterogeneity_evidence)
table(qc_summary$directional_pleiotropy_evidence)
table(qc_summary$presso_global_evidence)
#保存
write.csv(
  qc_summary,
  "05_results/06_summary/QC_summary.csv",
  row.names = FALSE
)
list.files("05_results/06_summary")
# ============================================================
# 16. Prepare forward IVW summary forest data
# ============================================================

forward_plot_dat <- ivw_summary %>%
  dplyr::filter(direction == "Forward") %>%
  dplyr::mutate(
    analysis_label = dplyr::case_when(
      analysis == "Primary GCST90018849" ~
        "Primary (GCST90018849)",
      
      analysis == "Threshold sensitivity P<=1e-5" ~
        "Threshold sensitivity (P <= 1e-5)",
      
      analysis == "FinnGen" ~
        "FinnGen",
      
      analysis == "East Asian GCST90018629" ~
        "East Asian (GCST90018629)",
      
      TRUE ~ analysis
    ),
    
    analysis_label = factor(
      analysis_label,
      levels = rev(c(
        "Primary (GCST90018849)",
        "Threshold sensitivity (P <= 1e-5)",
        "FinnGen",
        "East Asian (GCST90018629)"
      ))
    ),
    
    exposure_label = factor(
      exposure_label,
      levels = c(
        "Veillonella",
        "Veillonellaceae"
      )
    )
  )

forward_plot_dat %>%
  dplyr::select(
    exposure_label,
    analysis_label,
    nsnp,
    estimate,
    ci_lower,
    ci_upper,
    pval
  )
# ============================================================
# 17. Publication-level forward IVW forest plot
# ============================================================

forward_forest <- ggplot(
  forward_plot_dat,
  aes(
    x = estimate,
    y = analysis_label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  facet_wrap(
    ~ exposure_label,
    ncol = 1
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3)
  ) +
  labs(
    x = "Odds ratio (95% CI)",
    y = NULL,
    title = "Forward Mendelian randomization analyses"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0
    ),
    panel.grid.minor = element_blank()
  )

forward_forest
#保存
ggplot2::ggsave(
  "06_figures/06_summary/Forward_IVW_summary_forest.pdf",
  plot = forward_forest,
  width = 8,
  height = 7
)

ggplot2::ggsave(
  "06_figures/06_summary/Forward_IVW_summary_forest.png",
  plot = forward_forest,
  width = 8,
  height = 7,
  dpi = 600
)
# ============================================================
# 18. Reverse IVW summary forest
# ============================================================

reverse_plot_dat <- ivw_summary %>%
  dplyr::filter(direction == "Reverse") %>%
  dplyr::mutate(
    outcome_label = factor(
      outcome_label,
      levels = rev(c(
        "Veillonella",
        "Veillonellaceae"
      ))
    )
  )

reverse_forest <- ggplot(
  reverse_plot_dat,
  aes(
    x = estimate,
    y = outcome_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  labs(
    x = "Beta estimate (95% CI)",
    y = NULL,
    title = "Reverse Mendelian randomization analyses"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0
    ),
    panel.grid.minor = element_blank()
  )

reverse_forest
#保存
ggplot2::ggsave(
  "06_figures/06_summary/Reverse_IVW_summary_forest.pdf",
  plot = reverse_forest,
  width = 7,
  height = 4
)

ggplot2::ggsave(
  "06_figures/06_summary/Reverse_IVW_summary_forest.png",
  plot = reverse_forest,
  width = 7,
  height = 4,
  dpi = 600
)
#检查
list.files("06_figures/06_summary")

