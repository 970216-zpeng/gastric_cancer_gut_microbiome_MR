library(readr)
library(dplyr)
library(flextable)
library(officer)

# ============================================================
# 1. Read audited master data
# ============================================================

audit_path <- file.path(
  "07_tables",
  "01_audit_data",
  "Table2_core_MR_audit_master.csv"
)

stopifnot(file.exists(audit_path))

dat <- read_csv(audit_path, show_col_types = FALSE)

# Basic QC
stopifnot(nrow(dat) == 10)
stopifnot(all(dat$nsnp_check))

required_cols <- c(
  "analysis",
  "direction",
  "exposure",
  "outcome",
  "instrument_threshold",
  "nsnp_independent_count_from_singleSNP",
  "beta_output",
  "se_output",
  "p_recalculated_from_beta_se",
  "beta_ci95_low_recalculated",
  "beta_ci95_high_recalculated",
  "OR_recalculated_exp_beta",
  "OR_CI95_low_recalculated",
  "OR_CI95_high_recalculated"
)

stopifnot(all(required_cols %in% names(dat)))
# ============================================================
# 2. Build publication-layer Table 2
# ============================================================

tab <- dat %>%
  mutate(
    
    Analysis_display = case_when(
      analysis == "Primary" ~ "Primary GC",
      analysis == "Threshold sensitivity" ~ "Instrument-threshold sensitivity",
      analysis == "FinnGen sensitivity" ~ "Alternative European outcome",
      analysis == "East Asian sensitivity" ~ "Cross-ancestry outcome",
      analysis == "Reverse MR" ~ "Reverse MR",
      TRUE ~ analysis
    ),
    
    Outcome_analysis_display = case_when(
      analysis == "Primary" ~
        "GCST90018849 — Primary GC",
      
      analysis == "Threshold sensitivity" ~
        "GCST90018849 — Instrument-threshold sensitivity",
      
      analysis == "FinnGen sensitivity" ~
        "FinnGen gastric cancer — Alternative European outcome",
      
      analysis == "East Asian sensitivity" ~
        "GCST90018629 — Cross-ancestry outcome",
      
      analysis == "Reverse MR" ~
        outcome,
      
      TRUE ~ outcome
    ),
    
    Threshold_display = case_when(
      instrument_threshold == "P <= 5e-6" ~ "P ≤ 5 × 10⁻⁶",
      instrument_threshold == "P <= 1e-5" ~ "P ≤ 1 × 10⁻⁵",
      instrument_threshold == "P <= 5e-8" ~ "P ≤ 5 × 10⁻⁸",
      TRUE ~ instrument_threshold
    ),
    
    Effect_display = case_when(
      
      direction == "Forward" ~ sprintf(
        "OR %.2f (%.2f–%.2f)",
        OR_recalculated_exp_beta,
        OR_CI95_low_recalculated,
        OR_CI95_high_recalculated
      ),
      
      direction == "Reverse" ~ sprintf(
        "β %.4f (%.3f to %.3f)",
        beta_output,
        beta_ci95_low_recalculated,
        beta_ci95_high_recalculated
      )
    ),
    
    P_display = sprintf(
      "%.3f",
      p_recalculated_from_beta_se
    )
  ) %>%
  
  transmute(
    Direction = direction,
    Exposure = exposure,
    `Outcome / analysis` = Outcome_analysis_display,
    `Instrument threshold` = Threshold_display,
    SNPs = nsnp_independent_count_from_singleSNP,
    `IVW estimate (95% CI)` = Effect_display,
    `P value` = P_display
  )
View(tab)
# 2保存csv
dir.create(
  "07_tables/02_publication_data",
  recursive = TRUE, 
  showWarnings = FALSE
)

publication_path <- file.path(
  "07_tables",
  "02_publication_data",
  "Table2_publication.csv"
)

write_csv(tab, publication_path)

stopifnot(file.exists(publication_path))
#检查
publication_path
# ============================================================
# 4. Create manuscript Word table
# ============================================================

ft <- flextable(tab)

ft <- theme_booktabs(ft)

ft <- bold(
  ft,
  part = "header"
)

ft <- fontsize(
  ft,
  size = 9.5,
  part = "all"
)

ft <- align(
  ft,
  j = c("Direction", "SNPs", "P value"),
  align = "center",
  part = "all"
)

ft <- align(
  ft,
  j = c(
    "Instrument threshold",
    "IVW estimate (95% CI)"
  ),
  align = "center",
  part = "all"
)

ft <- valign(
  ft,
  valign = "center",
  part = "all"
)

ft <- autofit(ft)

ft <- set_caption(
  ft,
  caption = paste0(
    "Table 2. Summary of the core bidirectional ",
    "Mendelian randomization analyses"
  )
)

ft <- add_footer_lines(
  ft,
  values = c(
    paste0(
      "Abbreviations: CI, confidence interval; GC, gastric cancer; ",
      "IVW, inverse-variance weighted; MR, Mendelian randomization; ",
      "OR, odds ratio."
    ),
    paste0(
      "Forward MR estimates are presented as ORs for gastric cancer, ",
      "whereas reverse MR estimates are presented as β coefficients ",
      "for continuous microbiome outcomes."
    ),
    paste0(
      "The prespecified Bonferroni-adjusted significance threshold ",
      "for the two co-primary forward analyses was P < 0.025."
    ),
    paste0(
      "FinnGen was treated as an alternative European outcome ",
      "sensitivity analysis rather than an independent replication ",
      "because FinnGen contributed to GCST90018849."
    )
  )
)

dir.create(
  "07_tables/03_outputs",
  recursive = TRUE,
  showWarnings = FALSE
)

word_path <- file.path(
  "07_tables",
  "03_outputs",
  "Table2_core_bidirectional_MR.docx"
)

save_as_docx(
  "Table 2" = ft,
  path = word_path
)

stopifnot(file.exists(word_path))

word_path