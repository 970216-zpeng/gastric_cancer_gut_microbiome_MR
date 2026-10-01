# ============================================================
# Supplementary Figures S1-S6
# Independent reproducible pipeline
# ============================================================

# ============================================================
# 0. Packages
# ============================================================

library(here)
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(TwoSampleMR)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. Canonical analysis-source map
# ============================================================

source_map_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_analysis_source_map.csv"
)

stopifnot(
  file.exists(source_map_file)
)

fig_map <- read_csv(
  source_map_file,
  show_col_types = FALSE
)

stopifnot(
  nrow(fig_map) == 18,
  !anyDuplicated(fig_map$analysis_id)
)

# ============================================================
# 2. Locked final SNP counts
# ============================================================

expected_final_n <- c(
  
  core_primary_V = 4,
  core_primary_F = 9,
  
  core_p1e5_V = 8,
  core_p1e5_F = 19,
  
  finngen_V = 4,
  finngen_F = 9,
  
  eastasian_V = 4,
  eastasian_F = 5,
  
  reverse_V = 4,
  reverse_F = 4,
  
  swedish_primary_V = 4,
  swedish_primary_F = 3,
  
  dmp_primary_V = 4,
  dmp_primary_F = 6,
  
  swedish_p1e5_V = 10,
  swedish_p1e5_F = 10,
  
  dmp_p1e5_V = 6,
  dmp_p1e5_F = 12
)

stopifnot(
  setequal(
    names(expected_final_n),
    fig_map$analysis_id
  )
)

# ============================================================
# 3. Generic reader
# ============================================================

read_fig_file <- function(rel_path) {
  
  full_path <- here(rel_path)
  
  if (!file.exists(full_path)) {
    
    stop(
      "Missing source file: ",
      rel_path
    )
  }
  
  ext <- tolower(
    tools::file_ext(full_path)
  )
  
  if (ext == "csv") {
    
    x <- read_csv(
      full_path,
      show_col_types = FALSE
    )
    
  } else if (ext == "rds") {
    
    x <- readRDS(
      full_path
    )
    
  } else {
    
    stop(
      "Unsupported source type: ",
      full_path
    )
  }
  
  as.data.frame(x)
}

# ============================================================
# 4. Outcome formatter
# ============================================================

format_fig_outcome <- function(
    raw,
    outcome_name,
    outcome_id
) {
  
  raw <- as.data.frame(raw)
  
  already_formatted <- all(
    c(
      "SNP",
      "beta.outcome",
      "se.outcome",
      "effect_allele.outcome",
      "other_allele.outcome"
    ) %in% names(raw)
  )
  
  if (already_formatted) {
    
    if (!"id.outcome" %in% names(raw)) {
      raw$id.outcome <- outcome_id
    }
    
    if (!"outcome" %in% names(raw)) {
      raw$outcome <- outcome_name
    }
    
    if (!"mr_keep.outcome" %in% names(raw)) {
      raw$mr_keep.outcome <- TRUE
    }
    
    return(raw)
  }
  
  required_raw <- c(
    "SNP",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "pval"
  )
  
  missing_cols <- setdiff(
    required_raw,
    names(raw)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Outcome file missing: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  tibble(
    SNP =
      as.character(raw$SNP),
    
    effect_allele.outcome =
      as.character(raw$effect_allele),
    
    other_allele.outcome =
      as.character(raw$other_allele),
    
    beta.outcome =
      as.numeric(raw$beta),
    
    se.outcome =
      as.numeric(raw$se),
    
    pval.outcome =
      as.numeric(raw$pval),
    
    eaf.outcome =
      if(
        "eaf" %in% names(raw)
      ){
        as.numeric(raw$eaf)
      } else {
        NA_real_
      },
    
    outcome =
      outcome_name,
    
    id.outcome =
      outcome_id,
    
    mr_keep.outcome =
      TRUE
  )
}

# ============================================================
# 5. Independently load one final harmonised dataset
# ============================================================

fig_sets <- load_fig_analysis <- function(i) {
  
  meta <- fig_map[i, ]
  
  analysis_id <-
    as.character(
      meta$analysis_id
    )
  
  expected_n <-
    unname(
      expected_final_n[
        analysis_id
      ]
    )
  
  # ----------------------------------------------------------
  # Direct harmonised source
  # ----------------------------------------------------------
  
  if (
    meta$source_mode ==
    "direct_harmonised"
  ) {
    
    dat <- read_fig_file(
      meta$harmonised_path
    )
    
    # ----------------------------------------------------------
    # Re-harmonise frozen exposure/outcome files
    # ----------------------------------------------------------
    
  } else {
    
    exposure_dat <- read_fig_file(
      meta$exposure_path
    )
    
    outcome_raw <- read_fig_file(
      meta$outcome_path
    )
    
    if (!"id.exposure" %in%
        names(exposure_dat)) {
      
      exposure_dat$id.exposure <-
        analysis_id
    }
    
    if (!"exposure" %in%
        names(exposure_dat)) {
      
      exposure_dat$exposure <-
        as.character(
          meta$exposure
        )
    }
    
    outcome_dat <- format_fig_outcome(
      
      raw =
        outcome_raw,
      
      outcome_name =
        as.character(
          meta$outcome
        ),
      
      outcome_id =
        paste0(
          analysis_id,
          "_outcome"
        )
    )
    
    dat <- TwoSampleMR::harmonise_data(
      
      exposure_dat =
        exposure_dat,
      
      outcome_dat =
        outcome_dat,
      
      action = 3
    )
  }
  
  dat <- as.data.frame(dat)
  
  stopifnot(
    all(
      c(
        "SNP",
        "beta.exposure",
        "se.exposure",
        "beta.outcome",
        "se.outcome",
        "mr_keep"
      ) %in%
        names(dat)
    )
  )
  
  dat$mr_keep <-
    as.logical(
      dat$mr_keep
    )
  
  final <- dat %>%
    filter(
      mr_keep %in% TRUE
    )
  
  if (
    nrow(final) !=
    expected_n
  ) {
    
    stop(
      analysis_id,
      ": expected ",
      expected_n,
      " final SNPs; observed ",
      nrow(final)
    )
  }
  
  # Normalize metadata for plotting functions
  final$id.exposure <-
    analysis_id
  
  final$id.outcome <-
    paste0(
      analysis_id,
      "_outcome"
    )
  
  final$exposure <-
    as.character(
      meta$exposure
    )
  
  final$outcome <-
    as.character(
      meta$outcome
    )
  
  list(
    meta = meta,
    final = final
  )
}

# ============================================================
# 6. Load all analyses
# ============================================================

fig_sets <- map(
  seq_len(
    nrow(fig_map)
  ),
  load_fig_analysis
)

names(fig_sets) <-
  fig_map$analysis_id

# ============================================================
# 7. Global loading QC
# ============================================================

fig_load_audit <- map_dfr(
  names(fig_sets),
  function(id) {
    
    tibble(
      analysis_id = id,
      expected_n =
        expected_final_n[[id]],
      observed_n =
        nrow(
          fig_sets[[id]]$final
        )
    )
  }
)

stopifnot(
  nrow(fig_load_audit) == 18,
  all(
    fig_load_audit$expected_n ==
      fig_load_audit$observed_n
  ),
  sum(
    fig_load_audit$observed_n
  ) == 125
)

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY FIGURE DATA LOADING COMPLETE\n")
cat("====================================================\n")
cat("Analyses loaded: 18\n")
cat("Final SNP rows: 125\n")
cat("Final-count reconciliation: PASS\n")
cat("Independent clean loading: PASS\n")
cat("====================================================\n")
# ============================================================
# FIGURE S1
# Scatter plots of the primary forward MR analyses
# ============================================================

# ============================================================
# 8. Helper: run five MR estimators
# ============================================================

run_fig_mr <- function(dat) {
  
  set.seed(1234)
  
  TwoSampleMR::mr(
    dat,
    method_list = c(
      "mr_ivw",
      "mr_weighted_median",
      "mr_egger_regression",
      "mr_weighted_mode",
      "mr_simple_mode"
    )
  )
}

# ============================================================
# 9. Primary forward datasets
# ============================================================

figS1_V_dat <-
  fig_sets[["core_primary_V"]]$final

figS1_F_dat <-
  fig_sets[["core_primary_F"]]$final

stopifnot(
  nrow(figS1_V_dat) == 4,
  nrow(figS1_F_dat) == 9
)

# ============================================================
# 10. MR estimates
# ============================================================

figS1_V_mr <-
  run_fig_mr(
    figS1_V_dat
  )

figS1_F_mr <-
  run_fig_mr(
    figS1_F_dat
  )

# ============================================================
# 11. Generate scatter plots
# ============================================================

plot_V_raw <-
  TwoSampleMR::mr_scatter_plot(
    figS1_V_mr,
    figS1_V_dat
  )[[1]]

plot_F_raw <-
  TwoSampleMR::mr_scatter_plot(
    figS1_F_mr,
    figS1_F_dat
  )[[1]]

# ============================================================
# 12. Publication formatting
# ============================================================

plot_V <- plot_V_raw +
  
  labs(
    title =
      "a  Veillonella",
    x =
      "SNP effect on microbial abundance",
    y =
      "SNP effect on gastric cancer"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0
      ),
    
    legend.title =
      element_blank(),
    
    legend.position =
      "bottom"
  )


plot_F <- plot_F_raw +
  
  labs(
    title =
      "b  Veillonellaceae",
    x =
      "SNP effect on microbial abundance",
    y =
      "SNP effect on gastric cancer"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0
      ),
    
    legend.title =
      element_blank(),
    
    legend.position =
      "bottom"
  )

# ============================================================
# 13. Combine
# ============================================================

figS1 <- (
  plot_V |
    plot_F
) +
  
  plot_annotation(
    title =
      "Scatter plots of the primary forward Mendelian randomization analyses"
  ) &
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )

# ============================================================
# 14. Output paths
# ============================================================

figS1_pdf <- here(
  "06_figures",
  "10_final",
  "FigureS1_primary_forward_scatter.pdf"
)

figS1_png <- here(
  "06_figures",
  "10_final",
  "FigureS1_primary_forward_scatter.png"
)

dir.create(
  dirname(figS1_pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 15. Save
# ============================================================

ggsave(
  figS1_pdf,
  figS1,
  width = 11,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  figS1_png,
  figS1,
  width = 11,
  height = 5.8,
  units = "in",
  dpi = 600
)

# ============================================================
# 16. Figure S1 QC
# ============================================================

stopifnot(
  file.exists(figS1_pdf),
  file.exists(figS1_png)
)

cat("\n")
cat("====================================================\n")
cat("FIGURE S1 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")
cat("Panel A SNPs: 4\n")
cat("Panel B SNPs: 9\n")
cat("Primary forward source QC: PASS\n")
cat("PDF + PNG output: PASS\n")
cat("====================================================\n")
names(fig_sets)
names(fig_sets[["core_primary_v"]])
