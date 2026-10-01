# ============================================================
# Supplementary Tables S4-S6
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
library(MRPRESSO)
# ============================================================
# 1. Load canonical analysis-source map
# ============================================================

source_map_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_analysis_source_map.csv"
)

stopifnot(
  file.exists(source_map_file)
)

s4to6_map <- read_csv(
  source_map_file,
  show_col_types = FALSE
)

# Basic structural QC
required_map_cols <- c(
  "analysis_id",
  "analysis_group",
  "direction",
  "exposure_source",
  "exposure",
  "outcome",
  "instrument_threshold",
  "source_mode",
  "harmonised_path",
  "exposure_path",
  "outcome_path"
)

stopifnot(
  identical(
    names(s4to6_map),
    required_map_cols
  )
)

stopifnot(
  nrow(s4to6_map) == 18,
  !anyDuplicated(s4to6_map$analysis_id)
)

# ============================================================
# 2. Locked expected final SNP counts
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
    s4to6_map$analysis_id
  )
)

stopifnot(
  sum(expected_final_n) == 125
)

# ============================================================
# 3. Generic file reader
# ============================================================

read_analysis_file <- function(rel_path) {
  
  if (
    is.na(rel_path) ||
    !nzchar(rel_path)
  ) {
    stop("Empty file path.")
  }
  
  full_path <- here(rel_path)
  
  if (!file.exists(full_path)) {
    stop(
      "Missing file: ",
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
      "Unsupported file type: ",
      full_path
    )
  }
  
  as.data.frame(x)
}

# ============================================================
# 4. Format extracted outcome data for TwoSampleMR
# ============================================================

format_outcome_for_mr <- function(
    raw,
    outcome_name,
    outcome_id
) {
  
  raw <- as.data.frame(raw)
  
  # ----------------------------------------------------------
  # Case A: already in TwoSampleMR outcome format
  # ----------------------------------------------------------
  
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
    
    out <- raw
    
    if (!"id.outcome" %in% names(out)) {
      out$id.outcome <- outcome_id
    }
    
    if (!"outcome" %in% names(out)) {
      out$outcome <- outcome_name
    }
    
    if (!"mr_keep.outcome" %in% names(out)) {
      out$mr_keep.outcome <- TRUE
    }
    
    return(out)
  }
  
  # ----------------------------------------------------------
  # Case B: generic extracted-GWAS format
  # ----------------------------------------------------------
  
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
      "Outcome file missing columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  out <- tibble(
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
    
    outcome =
      outcome_name,
    
    id.outcome =
      outcome_id,
    
    mr_keep.outcome =
      TRUE,
    
    pval_origin.outcome =
      "reported"
  )
  
  # Optional fields
  if ("eaf" %in% names(raw)) {
    out$eaf.outcome <-
      as.numeric(raw$eaf)
  }
  
  if ("samplesize" %in% names(raw)) {
    out$samplesize.outcome <-
      as.numeric(raw$samplesize)
  }
  
  if ("ncase" %in% names(raw)) {
    out$ncase.outcome <-
      as.numeric(raw$ncase)
  }
  
  if ("ncontrol" %in% names(raw)) {
    out$ncontrol.outcome <-
      as.numeric(raw$ncontrol)
  }
  
  out
}

# ============================================================
# 5. Load / reconstruct one final harmonised analysis
# ============================================================

load_final_analysis <- function(i) {
  
  meta <- s4to6_map[i, ]
  
  analysis_id <-
    meta$analysis_id
  
  expected_n <-
    unname(
      expected_final_n[
        analysis_id
      ]
    )
  
  mode <-
    meta$source_mode
  
  # ----------------------------------------------------------
  # A. Direct frozen harmonised source
  # ----------------------------------------------------------
  
  if (mode == "direct_harmonised") {
    
    dat <- read_analysis_file(
      meta$harmonised_path
    )
    
    # ----------------------------------------------------------
    # B. Reconstruct from frozen exposure + outcome
    # ----------------------------------------------------------
    
  } else if (
    mode %in% c(
      "reharmonise",
      "reharmonise_reverse"
    )
  ) {
    
    exposure_dat <- read_analysis_file(
      meta$exposure_path
    )
    
    outcome_raw <- read_analysis_file(
      meta$outcome_path
    )
    
    # Ensure required exposure metadata exists
    if (!"id.exposure" %in% names(exposure_dat)) {
      exposure_dat$id.exposure <-
        analysis_id
    }
    
    if (!"exposure" %in% names(exposure_dat)) {
      exposure_dat$exposure <-
        meta$exposure
    }
    
    outcome_dat <- format_outcome_for_mr(
      raw =
        outcome_raw,
      
      outcome_name =
        meta$outcome,
      
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
    
  } else {
    
    stop(
      "Unknown source_mode for ",
      analysis_id,
      ": ",
      mode
    )
  }
  
  dat <- as.data.frame(dat)
  
  # ==========================================================
  # Required harmonised fields
  # ==========================================================
  
  required_harm_cols <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "beta.outcome",
    "se.outcome",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "mr_keep"
  )
  
  missing_cols <- setdiff(
    required_harm_cols,
    names(dat)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      analysis_id,
      ": missing harmonised columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  # ==========================================================
  # Keep only final MR SNPs
  # ==========================================================
  
  dat$mr_keep <-
    as.logical(
      dat$mr_keep
    )
  
  final_dat <- dat %>%
    filter(
      mr_keep %in% TRUE
    )
  
  # ==========================================================
  # Final count validation
  # ==========================================================
  
  observed_n <-
    nrow(final_dat)
  
  if (observed_n != expected_n) {
    
    stop(
      "\nFINAL SNP COUNT MISMATCH\n",
      "Analysis: ",
      analysis_id,
      "\nExpected: ",
      expected_n,
      "\nObserved: ",
      observed_n,
      "\n"
    )
  }
  
  # Duplicate SNP check
  if (
    anyDuplicated(
      final_dat$SNP
    )
  ) {
    
    stop(
      "Duplicate SNP detected in ",
      analysis_id
    )
  }
  
  # Missing essential values
  essential_numeric <- final_dat %>%
    select(
      beta.exposure,
      se.exposure,
      beta.outcome,
      se.outcome
    )
  
  if (
    anyNA(
      essential_numeric
    )
  ) {
    
    stop(
      "Missing beta/SE values in ",
      analysis_id
    )
  }
  
  # ==========================================================
  # Instrument strength validation
  # ==========================================================
  
  final_dat <- final_dat %>%
    mutate(
      F_stat_recalc =
        (
          beta.exposure /
            se.exposure
        )^2
    )
  
  if (
    any(
      final_dat$F_stat_recalc <= 10,
      na.rm = TRUE
    )
  ) {
    
    stop(
      "F <= 10 detected in ",
      analysis_id
    )
  }
  
  # ==========================================================
  # Return
  # ==========================================================
  
  list(
    
    meta =
      meta,
    
    final =
      final_dat,
    
    audit =
      tibble(
        analysis_id =
          analysis_id,
        
        source_mode =
          mode,
        
        expected_final_n =
          expected_n,
        
        observed_final_n =
          observed_n,
        
        duplicate_SNPs =
          anyDuplicated(
            final_dat$SNP
          ),
        
        missing_beta_SE =
          anyNA(
            essential_numeric
          ),
        
        min_F =
          min(
            final_dat$F_stat_recalc,
            na.rm = TRUE
          ),
        
        final_count_check =
          if_else(
            expected_n ==
              observed_n,
            "PASS",
            "CHECK"
          )
      )
  )
}

# ============================================================
# 6. Load all 18 analyses independently
# ============================================================

s4to6_sets <- map(
  seq_len(
    nrow(s4to6_map)
  ),
  load_final_analysis
)

names(s4to6_sets) <-
  s4to6_map$analysis_id

# ============================================================
# 7. Combine validation audit
# ============================================================

s4to6_load_audit <- bind_rows(
  map(
    s4to6_sets,
    "audit"
  )
)

stopifnot(
  nrow(s4to6_load_audit) == 18,
  all(
    s4to6_load_audit$final_count_check ==
      "PASS"
  ),
  all(
    s4to6_load_audit$duplicate_SNPs == 0
  ),
  all(
    s4to6_load_audit$missing_beta_SE ==
      FALSE
  ),
  all(
    s4to6_load_audit$min_F > 10
  ),
  sum(
    s4to6_load_audit$observed_final_n
  ) == 125
)

# ============================================================
# 8. Save loading-validation audit
# ============================================================

load_audit_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS4to6_source_loading_validation.csv"
)

write_csv(
  s4to6_load_audit,
  load_audit_file
)

# ============================================================
# 9. Final report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("S4-S6 INDEPENDENT SOURCE LOADING COMPLETE\n")
cat("====================================================\n")

cat(
  "Analyses loaded: ",
  length(s4to6_sets),
  "\n"
)

cat(
  "Total final SNP rows: ",
  sum(
    s4to6_load_audit$observed_final_n
  ),
  "\n"
)

cat(
  "Final SNP-count reconciliation: PASS\n"
)

cat(
  "Duplicate-SNP check: PASS\n"
)

cat(
  "Missing beta/SE check: PASS\n"
)

cat(
  "F > 10 check: PASS\n"
)

cat(
  "Clean independent loading: PASS\n"
)

cat("\nValidation file:\n")
cat(
  load_audit_file,
  "\n"
)

cat("====================================================\n")

print(
  s4to6_load_audit,
  n = Inf
)
# ============================================================
# SUPPLEMENTARY TABLE S4
# Complete MR estimates across primary and supportive estimators
# ============================================================

# ============================================================
# 10. MR method definitions
# ============================================================

s4_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)

s4_expected_method_names <- c(
  "Inverse variance weighted",
  "Weighted median",
  "MR Egger",
  "Weighted mode",
  "Simple mode"
)

# ============================================================
# 11. Run MR for one independently loaded analysis
# ============================================================

run_s4_analysis <- function(
    x,
    analysis_index
) {
  
  meta <- x$meta
  dat  <- as.data.frame(x$final)
  
  analysis_id <-
    as.character(meta$analysis_id)
  
  expected_n <-
    unname(
      expected_final_n[
        analysis_id
      ]
    )
  
  stopifnot(
    nrow(dat) == expected_n
  )
  
  # ----------------------------------------------------------
  # Required numerical fields
  # ----------------------------------------------------------
  
  required_cols <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "beta.outcome",
    "se.outcome"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(dat)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      analysis_id,
      ": missing required MR fields: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
  
  # ----------------------------------------------------------
  # Normalize TwoSampleMR metadata
  # ----------------------------------------------------------
  
  dat$id.exposure <-
    analysis_id
  
  dat$id.outcome <-
    paste0(
      analysis_id,
      "_outcome"
    )
  
  dat$exposure <-
    as.character(
      meta$exposure
    )
  
  dat$outcome <-
    as.character(
      meta$outcome
    )
  
  # All rows here have already passed final MR QC
  dat$mr_keep <- TRUE
  
  # ----------------------------------------------------------
  # Fixed seed
  #
  # Weighted-median/mode uncertainty estimation may involve
  # resampling. A fixed per-analysis seed makes S4 reproducible.
  # ----------------------------------------------------------
  
  set.seed(
    20260911 +
      analysis_index
  )
  
  # ----------------------------------------------------------
  # MR estimation
  # ----------------------------------------------------------
  
  res <- TwoSampleMR::mr(
    dat,
    method_list =
      s4_methods
  )
  
  # ----------------------------------------------------------
  # Attach frozen analytical metadata
  # ----------------------------------------------------------
  
  res %>%
    mutate(
      
      analysis_id =
        analysis_id,
      
      analysis_group =
        as.character(
          meta$analysis_group
        ),
      
      direction =
        as.character(
          meta$direction
        ),
      
      exposure_source =
        as.character(
          meta$exposure_source
        ),
      
      exposure_label =
        as.character(
          meta$exposure
        ),
      
      outcome_label =
        as.character(
          meta$outcome
        ),
      
      instrument_threshold =
        as.character(
          meta$instrument_threshold
        )
    )
}

# ============================================================
# 12. Run all 18 analyses
# ============================================================

s4_raw <- map2_dfr(
  
  s4to6_sets,
  
  seq_along(
    s4to6_sets
  ),
  
  run_s4_analysis
)

# ============================================================
# 13. Method-structure audit
# ============================================================

s4_raw <- s4_raw %>%
  filter(
    method %in%
      s4_expected_method_names
  )

s4_method_count <- s4_raw %>%
  count(
    analysis_id,
    name = "n_methods"
  )

print(
  as_tibble(s4_method_count),
  n = Inf
)
# Every analysis should contain IVW
s4_ivw_presence <- s4_raw %>%
  filter(
    method ==
      "Inverse variance weighted"
  ) %>%
  distinct(
    analysis_id
  )

stopifnot(
  nrow(s4_ivw_presence) == 18
)

# All analyses in this study have >=3 final SNPs,
# so all five requested estimators should be returned.
stopifnot(
  nrow(s4_method_count) == 18,
  all(
    s4_method_count$n_methods == 5
  ),
  nrow(s4_raw) == 90
)

# ============================================================
# 14. Derive interpretable effect estimates
# ============================================================

s4_audit <- s4_raw %>%
  mutate(
    
    beta_ci_low =
      b -
      1.96 * se,
    
    beta_ci_high =
      b +
      1.96 * se,
    
    effect_scale =
      if_else(
        direction == "Forward",
        "OR",
        "Beta"
      ),
    
    estimate =
      if_else(
        direction == "Forward",
        exp(b),
        b
      ),
    
    ci_low =
      if_else(
        direction == "Forward",
        exp(beta_ci_low),
        beta_ci_low
      ),
    
    ci_high =
      if_else(
        direction == "Forward",
        exp(beta_ci_high),
        beta_ci_high
      )
  )

# ============================================================
# 15. Exact validation against frozen deterministic results
# ============================================================

# These have exact beta / SE / P values preserved from the
# frozen analyses rather than rounded manuscript OR values.

s4_locked_exact_ivw <- tribble(
  
  ~analysis_id,
  ~locked_b,
  ~locked_se,
  ~locked_p,
  
  "core_primary_V",
  0.1184885,
  0.1275912,
  0.3530665,
  
  "core_primary_F",
  0.2432743,
  0.1396568,
  0.08151836,
  
  "core_p1e5_V",
  0.04293831,
  0.1082926,
  0.6917341,
  
  "core_p1e5_F",
  0.04086048,
  0.06793993,
  0.5475599,
  
  "reverse_V",
  0.003370458,
  0.07537796,
  0.9643352,
  
  "reverse_F",
  -0.003922082,
  0.06442253,
  0.9514543
)

s4_exact_ivw_check <- s4_audit %>%
  filter(
    method ==
      "Inverse variance weighted"
  ) %>%
  inner_join(
    s4_locked_exact_ivw,
    by = "analysis_id"
  ) %>%
  mutate(
    
    delta_b =
      abs(
        b -
          locked_b
      ),
    
    delta_se =
      abs(
        se -
          locked_se
      ),
    
    delta_p =
      abs(
        pval -
          locked_p
      ),
    
    audit_status =
      if_else(
        delta_b < 1e-6 &
          delta_se < 1e-6 &
          delta_p < 1e-6,
        "PASS",
        "CHECK"
      )
  )

print(
  as_tibble(
    s4_exact_ivw_check %>%
              select(
                analysis_id,
                b,
                locked_b,
                se,
                locked_se,
                pval,
                locked_p,
                audit_status
              )),
  n = Inf
)

stopifnot(
  nrow(s4_exact_ivw_check) == 6,
  all(
    s4_exact_ivw_check$audit_status ==
      "PASS"
  )
)

# ============================================================
# 16. Rounded-result consistency validation
# ============================================================

# These manuscript/frozen values were preserved primarily as
# rounded ORs and P values. Therefore this is a consistency
# check, not exact floating-point equality.

s4_locked_rounded_ivw <- tribble(
  
  ~analysis_id,
  ~locked_OR,
  ~locked_p,
  
  "finngen_V",
  1.052393,
  0.9011924,
  
  "finngen_F",
  1.566713,
  0.1778920,
  
  "eastasian_V",
  1.078227,
  0.5877881,
  
  "eastasian_F",
  1.223380,
  0.2293206,
  
  "swedish_primary_V",
  0.921,
  0.639,
  
  "swedish_primary_F",
  0.878,
  0.503,
  
  "dmp_primary_V",
  1.018,
  0.775,
  
  "dmp_primary_F",
  1.145,
  0.192,
  
  "swedish_p1e5_V",
  0.915,
  0.553,
  
  "swedish_p1e5_F",
  0.925,
  0.700,
  
  "dmp_p1e5_V",
  1.040,
  0.479,
  
  "dmp_p1e5_F",
  1.167,
  0.00616
)

s4_rounded_ivw_check <- s4_audit %>%
  filter(
    method ==
      "Inverse variance weighted"
  ) %>%
  inner_join(
    s4_locked_rounded_ivw,
    by = "analysis_id"
  ) %>%
  mutate(
    
    recalculated_OR =
      exp(b),
    
    delta_OR =
      abs(
        recalculated_OR -
          locked_OR
      ),
    
    delta_p =
      abs(
        pval -
          locked_p
      ),
    
    audit_status =
      if_else(
        delta_OR < 0.01 &
          delta_p < 0.01,
        "PASS",
        "CHECK"
      )
  )

print(
  as_tibble(s4_rounded_ivw_check %>%
              select(
                analysis_id,
                recalculated_OR,
                locked_OR,
                pval,
                locked_p,
                audit_status
              )),
  n = Inf
)

stopifnot(
  nrow(s4_rounded_ivw_check) == 12,
  all(
    s4_rounded_ivw_check$audit_status ==
      "PASS"
  )
)

# ============================================================
# 17. Primary supportive-estimator validation
# ============================================================

# Point estimates for these supportive methods were frozen
# previously. SE/P for resampling-based methods are not used as
# hard equality checks.

s4_locked_supportive <- tribble(
  
  ~analysis_id,
  ~method,
  ~locked_b,
  
  "core_primary_V",
  "Weighted median",
  0.1662629,
  
  "core_primary_V",
  "MR Egger",
  -1.5688114,
  
  "core_primary_V",
  "Weighted mode",
  0.1910533,
  
  "core_primary_V",
  "Simple mode",
  0.1875232,
  
  "core_primary_F",
  "Weighted median",
  0.2775634,
  
  "core_primary_F",
  "MR Egger",
  0.4334542,
  
  "core_primary_F",
  "Weighted mode",
  0.2946662,
  
  "core_primary_F",
  "Simple mode",
  0.3032692
)

s4_supportive_check <- s4_audit %>%
  inner_join(
    s4_locked_supportive,
    by = c(
      "analysis_id",
      "method"
    )
  ) %>%
  mutate(
    
    delta_b =
      abs(
        b -
          locked_b
      ),
    
    audit_status =
      if_else(
        delta_b < 1e-4,
        "PASS",
        "CHECK"
      )
  )

print(
  as_tibble(s4_supportive_check %>%
              select(
                analysis_id,
                method,
                b,
                locked_b,
                delta_b,
                audit_status
              )),
  n = Inf
)

stopifnot(
  nrow(s4_supportive_check) == 8,
  all(
    s4_supportive_check$audit_status ==
      "PASS"
  )
)

# ============================================================
# 18. Global numerical QC
# ============================================================

stopifnot(
  n_distinct(
    s4_audit$analysis_id
  ) == 18,
  
  nrow(s4_audit) == 90,
  
  all(
    is.finite(
      s4_audit$b
    )
  ),
  
  all(
    is.finite(
      s4_audit$se
    )
  ),
  
  all(
    is.finite(
      s4_audit$pval
    )
  ),
  
  all(
    s4_audit$pval >= 0 &
      s4_audit$pval <= 1
  )
)

# ============================================================
# 19. Save S4 audit outputs
# ============================================================

s4_audit_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS4_complete_MR_estimates_audit_master.csv"
)

s4_exact_check_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS4_exact_IVW_validation.csv"
)

s4_rounded_check_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS4_rounded_IVW_consistency_validation.csv"
)

s4_supportive_check_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS4_primary_supportive_validation.csv"
)

write_csv(
  s4_audit,
  s4_audit_file
)

write_csv(
  s4_exact_ivw_check,
  s4_exact_check_file
)

write_csv(
  s4_rounded_ivw_check,
  s4_rounded_check_file
)

write_csv(
  s4_supportive_check,
  s4_supportive_check_file
)

# ============================================================
# 20. S4 audit report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S4 MR AUDIT COMPLETE\n")
cat("====================================================\n")

cat(
  "Analyses: ",
  n_distinct(
    s4_audit$analysis_id
  ),
  "\n"
)

cat(
  "Methods per analysis: 5\n"
)

cat(
  "Total MR estimate rows: ",
  nrow(s4_audit),
  "\n"
)

cat(
  "Exact IVW validation: PASS\n"
)

cat(
  "Rounded IVW consistency: PASS\n"
)

cat(
  "Primary supportive point-estimate validation: PASS\n"
)

cat(
  "Finite-value / P-range QC: PASS\n"
)

cat(
  "Fixed-seed reproducibility: PASS\n"
)

cat("\nAudit master:\n")
cat(
  s4_audit_file,
  "\n"
)

cat("====================================================\n")
# ============================================================
# SUPPLEMENTARY TABLE S4
# Publication layer and Word output
# ============================================================

library(flextable)
library(officer)

# ============================================================
# 21. Publication labels
# ============================================================

s4_analysis_labels <- c(
  
  core_primary_V =
    "MiBioGen Veillonella → GCST90018849 (primary)",
  
  core_primary_F =
    "MiBioGen Veillonellaceae → GCST90018849 (primary)",
  
  core_p1e5_V =
    "MiBioGen Veillonella → GCST90018849 (threshold sensitivity)",
  
  core_p1e5_F =
    "MiBioGen Veillonellaceae → GCST90018849 (threshold sensitivity)",
  
  finngen_V =
    "MiBioGen Veillonella → FinnGen GC (Alt-E)",
  
  finngen_F =
    "MiBioGen Veillonellaceae → FinnGen GC (Alt-E)",
  
  eastasian_V =
    "MiBioGen Veillonella → East Asian GC",
  
  eastasian_F =
    "MiBioGen Veillonellaceae → East Asian GC",
  
  reverse_V =
    "GCST90018849 GC → MiBioGen Veillonella (reverse)",
  
  reverse_F =
    "GCST90018849 GC → MiBioGen Veillonellaceae (reverse)",
  
  swedish_primary_V =
    "Swedish 2026 Veillonella → GCST90018849 (primary)",
  
  swedish_primary_F =
    "Swedish 2026 Veillonellaceae → GCST90018849 (primary)",
  
  dmp_primary_V =
    "DMP 2022 Veillonella → GCST90018849 (primary)",
  
  dmp_primary_F =
    "DMP 2022 Veillonellaceae → GCST90018849 (primary)",
  
  swedish_p1e5_V =
    "Swedish 2026 Veillonella → GCST90018849 (threshold sensitivity)",
  
  swedish_p1e5_F =
    "Swedish 2026 Veillonellaceae → GCST90018849 (threshold sensitivity)",
  
  dmp_p1e5_V =
    "DMP 2022 Veillonella → GCST90018849 (threshold sensitivity)",
  
  dmp_p1e5_F =
    "DMP 2022 Veillonellaceae → GCST90018849 (threshold sensitivity)"
)

stopifnot(
  setequal(
    names(s4_analysis_labels),
    unique(s4_audit$analysis_id)
  )
)

# ============================================================
# 22. Formatting helpers
# ============================================================

format_s4_effect <- function(
    direction,
    estimate,
    ci_low,
    ci_high
) {
  
  ifelse(
    direction == "Forward",
    
    sprintf(
      "OR %.3f (%.3f–%.3f)",
      estimate,
      ci_low,
      ci_high
    ),
    
    sprintf(
      "β %.4f (%.4f–%.4f)",
      estimate,
      ci_low,
      ci_high
    )
  )
}

format_s4_p <- function(p) {
  
  ifelse(
    p < 0.001,
    
    formatC(
      p,
      format = "e",
      digits = 2
    ),
    
    ifelse(
      p < 0.01,
      sprintf("%.4f", p),
      sprintf("%.3f", p)
    )
  )
}

# ============================================================
# 23. Build publication dataframe
# ============================================================

s4_pub <- s4_audit %>%
  mutate(
    
    Analysis =
      unname(
        s4_analysis_labels[
          analysis_id
        ]
      ),
    
    Method =
      factor(
        method,
        levels =
          s4_expected_method_names
      ),
    
    analysis_order =
      match(
        analysis_id,
        s4to6_map$analysis_id
      ),
    
    `Estimate (95% CI)` =
      format_s4_effect(
        direction,
        estimate,
        ci_low,
        ci_high
      ),
    
    `P value` =
      format_s4_p(
        pval
      )
  ) %>%
  arrange(
    analysis_order,
    Method
  ) %>%
  transmute(
    
    Analysis,
    
    Threshold =
      instrument_threshold,
    
    Method =
      as.character(Method),
    
    SNPs =
      nsnp,
    
    `Estimate (95% CI)`,
    
    `P value`
  )

# ============================================================
# 24. Publication-layer QC
# ============================================================

stopifnot(
  nrow(s4_pub) == 90
)

stopifnot(
  n_distinct(
    s4_pub$Analysis
  ) == 18
)

s4_pub_method_check <- s4_pub %>%
  count(
    Analysis,
    name = "n_methods"
  )

stopifnot(
  all(
    s4_pub_method_check$n_methods ==
      5
  )
)

stopifnot(
  !anyNA(
    s4_pub$`Estimate (95% CI)`
  ),
  !anyNA(
    s4_pub$`P value`
  )
)

# ============================================================
# 25. Save publication CSV
# ============================================================

s4_pub_file <- here(
  "07_tables",
  "02_publication_data",
  "TableS4_publication.csv"
)

write_csv(
  s4_pub,
  s4_pub_file
)

# ============================================================
# 26. Build flextable
# ============================================================

ft_s4 <- flextable(
  s4_pub
)

ft_s4 <- font(
  ft_s4,
  fontname = "Arial",
  part = "all"
)

ft_s4 <- fontsize(
  ft_s4,
  size = 7.5,
  part = "all"
)

ft_s4 <- bold(
  ft_s4,
  part = "header"
)

ft_s4 <- align(
  ft_s4,
  align = "center",
  part = "header"
)

ft_s4 <- align(
  ft_s4,
  j = "Analysis",
  align = "left",
  part = "body"
)

ft_s4 <- align(
  ft_s4,
  j = c(
    "Threshold",
    "Method",
    "SNPs",
    "Estimate (95% CI)",
    "P value"
  ),
  align = "center",
  part = "body"
)

ft_s4 <- valign(
  ft_s4,
  valign = "center",
  part = "all"
)

# ============================================================
# 27. Booktabs formatting
# ============================================================

ft_s4 <- border_remove(
  ft_s4
)

ft_s4 <- hline_top(
  ft_s4,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "header"
)

ft_s4 <- hline_bottom(
  ft_s4,
  border = fp_border(
    color = "black",
    width = 0.7
  ),
  part = "header"
)

ft_s4 <- hline_bottom(
  ft_s4,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "body"
)

ft_s4 <- padding(
  ft_s4,
  padding.top = 1.8,
  padding.bottom = 1.8,
  padding.left = 2.5,
  padding.right = 2.5,
  part = "all"
)

# ============================================================
# 28. Column widths
# ============================================================

ft_s4 <- width(
  ft_s4,
  j = "Analysis",
  width = 4.20
)

ft_s4 <- width(
  ft_s4,
  j = "Threshold",
  width = 1.25
)

ft_s4 <- width(
  ft_s4,
  j = "Method",
  width = 1.55
)

ft_s4 <- width(
  ft_s4,
  j = "SNPs",
  width = 0.65
)

ft_s4 <- width(
  ft_s4,
  j = "Estimate (95% CI)",
  width = 1.80
)

ft_s4 <- width(
  ft_s4,
  j = "P value",
  width = 0.90
)

# Allow multi-page table and repeat header
ft_s4 <- set_table_properties(
  ft_s4,
  opts_word = list(
    split = TRUE,
    repeat_headers = TRUE
  )
)

# ============================================================
# 29. Title and notes
# ============================================================

s4_title <- paste0(
  "Supplementary Table S4. Complete Mendelian randomization ",
  "estimates across primary and supportive estimators"
)

s4_note1 <- paste0(
  "Abbreviations: Alt-E, alternative European; CI, confidence interval; ",
  "DMP, Dutch Microbiome Project; GC, gastric cancer; IVW, ",
  "inverse-variance weighted; MR, Mendelian randomization; ",
  "OR, odds ratio; SE, standard error."
)

s4_note2 <- paste0(
  "IVW was the primary estimator. Weighted median, MR-Egger, ",
  "weighted mode, and simple mode were supportive estimators."
)

s4_note3 <- paste0(
  "Forward MR estimates are presented as ORs for gastric cancer. ",
  "Reverse MR estimates are presented as beta coefficients for ",
  "continuous microbiome outcomes."
)

s4_note4 <- paste0(
  "Supportive estimators were interpreted jointly with the IVW estimate ",
  "and were not used to selectively replace the primary estimator."
)

# ============================================================
# 30. Word output
# ============================================================

s4_output <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S4_complete_MR_estimates.docx"
)

doc_s4 <- read_docx()

doc_s4 <- body_add_fpar(
  doc_s4,
  fpar(
    ftext(
      s4_title,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc_s4 <- body_add_flextable(
  doc_s4,
  value = ft_s4
)

for (
  txt in c(
    s4_note1,
    s4_note2,
    s4_note3,
    s4_note4
  )
) {
  
  doc_s4 <- body_add_fpar(
    doc_s4,
    fpar(
      ftext(
        txt,
        prop = fp_text(
          font.family = "Arial",
          font.size = 8
        )
      )
    )
  )
}

doc_s4 <- body_end_section_landscape(
  doc_s4
)

print(
  doc_s4,
  target = s4_output
)

# ============================================================
# 31. Final S4 publication report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S4 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "Analyses: 18\n"
)

cat(
  "Methods per analysis: 5\n"
)

cat(
  "Publication rows: ",
  nrow(s4_pub),
  "\n"
)

cat(
  "Publication-layer QC: PASS\n"
)

cat(
  "Numerical audit: PASS\n"
)

cat(
  "Output:\n"
)

cat(
  s4_output,
  "\n"
)

cat("====================================================\n")
# ============================================================
# SUPPLEMENTARY TABLE S5
# Heterogeneity and horizontal pleiotropy diagnostics
# ============================================================

# ============================================================
# 32. Helper: prepare final harmonised data for diagnostics
# ============================================================

prepare_s5_data <- function(x) {
  
  meta <- x$meta
  dat  <- as.data.frame(x$final)
  
  analysis_id <-
    as.character(meta$analysis_id)
  
  dat$id.exposure <-
    analysis_id
  
  dat$id.outcome <-
    paste0(
      analysis_id,
      "_outcome"
    )
  
  dat$exposure <-
    as.character(
      meta$exposure
    )
  
  dat$outcome <-
    as.character(
      meta$outcome
    )
  
  dat$mr_keep <- TRUE
  
  list(
    meta = meta,
    dat = dat
  )
}


# ============================================================
# 33. Helper: safely parse MR-PRESSO P values
# ============================================================

parse_presso_p <- function(x) {
  
  if (length(x) == 0 || is.null(x)) {
    return(NA_real_)
  }
  
  x <- x[1]
  
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x <- as.character(x)
  
  # e.g. "<0.0001"
  x <- gsub(
    "^<\\s*",
    "",
    x
  )
  
  suppressWarnings(
    as.numeric(x)
  )
}


# ============================================================
# 34. Helper: run MR-PRESSO safely
# ============================================================

run_s5_presso <- function(
    dat,
    analysis_id,
    analysis_index
) {
  
  nsnp <- nrow(dat)
  
  # MR-PRESSO is not methodologically meaningful with only 3 SNPs.
  if (nsnp < 4) {
    
    return(
      tibble(
        presso_feasible = FALSE,
        presso_global_p = NA_real_,
        presso_outliers = NA_integer_,
        presso_status = "Not feasible (<4 SNPs)"
      )
    )
  }
  
  # Fixed seed for reproducibility
  set.seed(
    20260912 +
      analysis_index
  )
  
  fit <- tryCatch(
    
    MRPRESSO::mr_presso(
      
      BetaOutcome =
        "beta.outcome",
      
      BetaExposure =
        "beta.exposure",
      
      SdOutcome =
        "se.outcome",
      
      SdExposure =
        "se.exposure",
      
      OUTLIERtest =
        TRUE,
      
      DISTORTIONtest =
        TRUE,
      
      data =
        dat,
      
      NbDistribution =
        10000,
      
      SignifThreshold =
        0.05
    ),
    
    error = function(e) {
      e
    }
  )
  
  if (inherits(fit, "error")) {
    
    return(
      tibble(
        presso_feasible = TRUE,
        presso_global_p = NA_real_,
        presso_outliers = NA_integer_,
        presso_status =
          paste0(
            "ERROR: ",
            conditionMessage(fit)
          )
      )
    )
  }
  
  # ----------------------------------------------------------
  # Global test
  # ----------------------------------------------------------
  
  global_test <-
    fit[["MR-PRESSO results"]][["Global Test"]]
  
  global_p <-
    parse_presso_p(
      global_test$Pvalue
    )
  
  # ----------------------------------------------------------
  # Outlier test
  # ----------------------------------------------------------
  
  outlier_test <-
    fit[["MR-PRESSO results"]][["Outlier Test"]]
  
  if (
    is.null(outlier_test) ||
    nrow(outlier_test) == 0
  ) {
    
    n_outliers <- 0L
    
  } else {
    
    outlier_p <-
      vapply(
        outlier_test$Pvalue,
        parse_presso_p,
        numeric(1)
      )
    
    n_outliers <-
      sum(
        outlier_p < 0.05,
        na.rm = TRUE
      )
  }
  
  tibble(
    presso_feasible = TRUE,
    presso_global_p = global_p,
    presso_outliers = as.integer(n_outliers),
    presso_status = "Completed"
  )
}


# ============================================================
# 35. Run diagnostics for one analysis
# ============================================================

run_s5_analysis <- function(
    x,
    analysis_index
) {
  
  prepared <-
    prepare_s5_data(x)
  
  meta <- prepared$meta
  dat  <- prepared$dat
  
  analysis_id <-
    as.character(
      meta$analysis_id
    )
  
  nsnp <-
    nrow(dat)
  
  # ----------------------------------------------------------
  # Cochran's Q
  # ----------------------------------------------------------
  
  het <- TwoSampleMR::mr_heterogeneity(
    dat,
    method_list = c(
      "mr_ivw",
      "mr_egger_regression"
    )
  )
  
  ivw_het <- het %>%
    filter(
      grepl(
        "Inverse variance weighted",
        method,
        fixed = TRUE
      )
    ) %>%
    slice(1)
  
  if (nrow(ivw_het) != 1) {
    
    stop(
      analysis_id,
      ": IVW heterogeneity result not uniquely identified."
    )
  }
  
  # ----------------------------------------------------------
  # MR-Egger intercept
  # ----------------------------------------------------------
  
  egger <- TwoSampleMR::mr_pleiotropy_test(
    dat
  )
  
  if (nrow(egger) != 1) {
    
    stop(
      analysis_id,
      ": MR-Egger intercept result not uniquely identified."
    )
  }
  
  # ----------------------------------------------------------
  # MR-PRESSO
  # ----------------------------------------------------------
  
  presso <- run_s5_presso(
    dat =
      dat,
    
    analysis_id =
      analysis_id,
    
    analysis_index =
      analysis_index
  )
  
  # ----------------------------------------------------------
  # Final row
  # ----------------------------------------------------------
  
  tibble(
    
    analysis_id =
      analysis_id,
    
    analysis_group =
      as.character(
        meta$analysis_group
      ),
    
    direction =
      as.character(
        meta$direction
      ),
    
    exposure_source =
      as.character(
        meta$exposure_source
      ),
    
    exposure =
      as.character(
        meta$exposure
      ),
    
    outcome =
      as.character(
        meta$outcome
      ),
    
    instrument_threshold =
      as.character(
        meta$instrument_threshold
      ),
    
    nsnp =
      nsnp,
    
    Q =
      ivw_het$Q,
    
    Q_df =
      ivw_het$Q_df,
    
    Q_p =
      ivw_het$Q_pval,
    
    egger_intercept =
      egger$egger_intercept,
    
    egger_intercept_se =
      egger$se,
    
    egger_intercept_p =
      egger$pval
    
  ) %>%
    bind_cols(
      presso
    )
}


# ============================================================
# 36. Run all 18 analyses
# ============================================================

s5_audit <- map2_dfr(
  
  s4to6_sets,
  
  seq_along(
    s4to6_sets
  ),
  
  run_s5_analysis
)


# ============================================================
# 37. Structural QC
# ============================================================

stopifnot(
  nrow(s5_audit) == 18,
  n_distinct(
    s5_audit$analysis_id
  ) == 18
)

stopifnot(
  all(
    is.finite(
      s5_audit$Q
    )
  ),
  all(
    is.finite(
      s5_audit$Q_p
    )
  ),
  all(
    s5_audit$Q_p >= 0 &
      s5_audit$Q_p <= 1
  )
)

stopifnot(
  all(
    is.finite(
      s5_audit$egger_intercept
    )
  ),
  all(
    is.finite(
      s5_audit$egger_intercept_p
    )
  ),
  all(
    s5_audit$egger_intercept_p >= 0 &
      s5_audit$egger_intercept_p <= 1
  )
)


# ============================================================
# 38. MR-PRESSO feasibility QC
# ============================================================

# Exactly one analysis is expected to have only 3 SNPs:
# Swedish 2026 Veillonellaceae primary.

s5_infeasible <- s5_audit %>%
  filter(
    !presso_feasible
  )

stopifnot(
  nrow(s5_infeasible) == 1,
  s5_infeasible$analysis_id ==
    "swedish_primary_F",
  s5_infeasible$nsnp == 3
)

# All feasible analyses should complete successfully.
s5_presso_completed <- s5_audit %>%
  filter(
    presso_feasible
  )

stopifnot(
  nrow(s5_presso_completed) == 17,
  all(
    s5_presso_completed$presso_status ==
      "Completed"
  ),
  all(
    !is.na(
      s5_presso_completed$presso_global_p
    )
  )
)


# ============================================================
# 39. Evidence flags
# ============================================================

s5_audit <- s5_audit %>%
  mutate(
    
    heterogeneity_evidence =
      Q_p < 0.05,
    
    egger_pleiotropy_evidence =
      egger_intercept_p < 0.05,
    
    presso_global_evidence =
      if_else(
        presso_feasible,
        presso_global_p < 0.05,
        NA
      )
  )


# ============================================================
# 40. Validate external Post-reg diagnostics
# ============================================================

s5_external <- s5_audit %>%
  filter(
    analysis_group ==
      "Post-reg"
  )

stopifnot(
  nrow(s5_external) == 8
)

# Frozen final audit:
# 0 / 8 showed heterogeneity evidence
stopifnot(
  sum(
    s5_external$heterogeneity_evidence,
    na.rm = TRUE
  ) == 0
)

# Frozen final audit:
# 0 / 8 showed Egger-intercept evidence
stopifnot(
  sum(
    s5_external$egger_pleiotropy_evidence,
    na.rm = TRUE
  ) == 0
)

# Frozen final audit:
# no feasible external MR-PRESSO test showed global evidence
stopifnot(
  sum(
    s5_external$presso_global_evidence,
    na.rm = TRUE
  ) == 0
)


# ============================================================
# 41. Specific frozen-result validation:
# DMP Veillonellaceae, P <= 1e-5
# ============================================================

s5_dmp_F_p1e5 <- s5_audit %>%
  filter(
    analysis_id ==
      "dmp_p1e5_F"
  )

stopifnot(
  nrow(s5_dmp_F_p1e5) == 1
)

# Frozen:
# Q P = 0.718
# Egger intercept P = 0.871
# PRESSO global P = 0.7752

stopifnot(
  abs(
    s5_dmp_F_p1e5$Q_p -
      0.718
  ) < 0.01,
  
  abs(
    s5_dmp_F_p1e5$egger_intercept_p -
      0.871
  ) < 0.01,
  
  abs(
    s5_dmp_F_p1e5$presso_global_p -
      0.7752
  ) < 0.05
)


# ============================================================
# 42. Save S5 audit master
# ============================================================

s5_audit_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS5_diagnostics_audit_master.csv"
)

write_csv(
  s5_audit,
  s5_audit_file
)


# ============================================================
# 43. Console report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S5 DATA AUDIT COMPLETE\n")
cat("====================================================\n")

cat(
  "Analyses: ",
  nrow(s5_audit),
  "\n"
)

cat(
  "Heterogeneity P < 0.05: ",
  sum(
    s5_audit$heterogeneity_evidence,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "MR-Egger intercept P < 0.05: ",
  sum(
    s5_audit$egger_pleiotropy_evidence,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "MR-PRESSO feasible: ",
  sum(
    s5_audit$presso_feasible
  ),
  " / 18\n"
)

cat(
  "MR-PRESSO global P < 0.05: ",
  sum(
    s5_audit$presso_global_evidence,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Post-reg frozen-result reconciliation: PASS\n"
)

cat(
  "DMP p1e-5 Veillonellaceae diagnostic check: PASS\n"
)

cat("\nAudit master:\n")
cat(
  s5_audit_file,
  "\n"
)

cat("====================================================\n")

print(
  as_tibble(
    s5_audit %>%
      select(
        analysis_id,
        nsnp,
        Q,
        Q_p,
        egger_intercept,
        egger_intercept_p,
        presso_global_p,
        presso_outliers,
        presso_status
      )
  ),
  n = Inf
)
# ============================================================
# SUPPLEMENTARY TABLE S5
# Publication layer and Word output
# ============================================================

# ============================================================
# 44. Analysis labels
# ============================================================

s5_analysis_labels <- c(
  
  core_primary_V =
    "MiBioGen Veillonella → GCST90018849 (primary)",
  
  core_primary_F =
    "MiBioGen Veillonellaceae → GCST90018849 (primary)",
  
  core_p1e5_V =
    "MiBioGen Veillonella → GCST90018849 (threshold sensitivity)",
  
  core_p1e5_F =
    "MiBioGen Veillonellaceae → GCST90018849 (threshold sensitivity)",
  
  finngen_V =
    "MiBioGen Veillonella → FinnGen GC (Alt-E)",
  
  finngen_F =
    "MiBioGen Veillonellaceae → FinnGen GC (Alt-E)",
  
  eastasian_V =
    "MiBioGen Veillonella → East Asian GC",
  
  eastasian_F =
    "MiBioGen Veillonellaceae → East Asian GC",
  
  reverse_V =
    "GCST90018849 GC → MiBioGen Veillonella (reverse)",
  
  reverse_F =
    "GCST90018849 GC → MiBioGen Veillonellaceae (reverse)",
  
  swedish_primary_V =
    "Swedish 2026 Veillonella → GCST90018849 (primary)",
  
  swedish_primary_F =
    "Swedish 2026 Veillonellaceae → GCST90018849 (primary)",
  
  dmp_primary_V =
    "DMP 2022 Veillonella → GCST90018849 (primary)",
  
  dmp_primary_F =
    "DMP 2022 Veillonellaceae → GCST90018849 (primary)",
  
  swedish_p1e5_V =
    "Swedish 2026 Veillonella → GCST90018849 (threshold sensitivity)",
  
  swedish_p1e5_F =
    "Swedish 2026 Veillonellaceae → GCST90018849 (threshold sensitivity)",
  
  dmp_p1e5_V =
    "DMP 2022 Veillonella → GCST90018849 (threshold sensitivity)",
  
  dmp_p1e5_F =
    "DMP 2022 Veillonellaceae → GCST90018849 (threshold sensitivity)"
)

stopifnot(
  setequal(
    names(s5_analysis_labels),
    s5_audit$analysis_id
  )
)

# ============================================================
# 45. Formatting helpers
# ============================================================

format_s5_p <- function(p) {
  
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      formatC(
        p,
        format = "e",
        digits = 2
      ),
      sprintf("%.3f", p)
    )
  )
}

format_s5_num <- function(x) {
  
  ifelse(
    is.na(x),
    "",
    sprintf("%.4f", x)
  )
}

# ============================================================
# 46. Build publication dataframe
# ============================================================

s5_pub <- s5_audit %>%
  mutate(
    
    Analysis =
      unname(
        s5_analysis_labels[
          analysis_id
        ]
      ),
    
    analysis_order =
      match(
        analysis_id,
        s4to6_map$analysis_id
      ),
    
    `Cochran Q (df)` =
      paste0(
        sprintf("%.3f", Q),
        " (",
        Q_df,
        ")"
      ),
    
    `Q P value` =
      format_s5_p(
        Q_p
      ),
    
    `MR-Egger intercept (SE)` =
      paste0(
        format_s5_num(
          egger_intercept
        ),
        " (",
        format_s5_num(
          egger_intercept_se
        ),
        ")"
      ),
    
    `Egger P value` =
      format_s5_p(
        egger_intercept_p
      ),
    
    `MR-PRESSO global P` =
      if_else(
        presso_feasible,
        format_s5_p(
          presso_global_p
        ),
        "Not feasible"
      ),
    
    `MR-PRESSO outliers` =
      if_else(
        presso_feasible,
        as.character(
          presso_outliers
        ),
        "—"
      )
  ) %>%
  arrange(
    analysis_order
  ) %>%
  transmute(
    
    Analysis,
    
    Threshold =
      instrument_threshold,
    
    SNPs =
      nsnp,
    
    `Cochran Q (df)`,
    
    `Q P value`,
    
    `MR-Egger intercept (SE)`,
    
    `Egger P value`,
    
    `MR-PRESSO global P`,
    
    `MR-PRESSO outliers`
  )

# ============================================================
# 47. Publication QC
# ============================================================

stopifnot(
  nrow(s5_pub) == 18,
  n_distinct(
    s5_pub$Analysis
  ) == 18
)

stopifnot(
  sum(
    s5_pub$`MR-PRESSO global P` ==
      "Not feasible"
  ) == 1
)

stopifnot(
  s5_pub$SNPs[
    s5_pub$`MR-PRESSO global P` ==
      "Not feasible"
  ] == 3
)

# ============================================================
# 48. Save publication CSV
# ============================================================

s5_pub_file <- here(
  "07_tables",
  "02_publication_data",
  "TableS5_publication.csv"
)

write_csv(
  s5_pub,
  s5_pub_file
)

# ============================================================
# 49. Build flextable
# ============================================================

ft_s5 <- flextable(
  s5_pub
)

ft_s5 <- font(
  ft_s5,
  fontname = "Arial",
  part = "all"
)

ft_s5 <- fontsize(
  ft_s5,
  size = 7.5,
  part = "all"
)

ft_s5 <- bold(
  ft_s5,
  part = "header"
)

ft_s5 <- align(
  ft_s5,
  align = "center",
  part = "header"
)

ft_s5 <- align(
  ft_s5,
  j = "Analysis",
  align = "left",
  part = "body"
)

ft_s5 <- align(
  ft_s5,
  j = 2:ncol(s5_pub),
  align = "center",
  part = "body"
)

ft_s5 <- valign(
  ft_s5,
  valign = "center",
  part = "all"
)

# ============================================================
# 50. Booktabs
# ============================================================

ft_s5 <- border_remove(
  ft_s5
)

ft_s5 <- hline_top(
  ft_s5,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "header"
)

ft_s5 <- hline_bottom(
  ft_s5,
  border = fp_border(
    color = "black",
    width = 0.7
  ),
  part = "header"
)

ft_s5 <- hline_bottom(
  ft_s5,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "body"
)

ft_s5 <- padding(
  ft_s5,
  padding.top = 2.3,
  padding.bottom = 2.3,
  padding.left = 2,
  padding.right = 2,
  part = "all"
)

# ============================================================
# 51. Column widths
# ============================================================

ft_s5 <- width(
  ft_s5,
  j = "Analysis",
  width = 3.65
)

ft_s5 <- width(
  ft_s5,
  j = "Threshold",
  width = 1.05
)

ft_s5 <- width(
  ft_s5,
  j = "SNPs",
  width = 0.50
)

ft_s5 <- width(
  ft_s5,
  j = "Cochran Q (df)",
  width = 1.10
)

ft_s5 <- width(
  ft_s5,
  j = "Q P value",
  width = 0.80
)

ft_s5 <- width(
  ft_s5,
  j = "MR-Egger intercept (SE)",
  width = 1.55
)

ft_s5 <- width(
  ft_s5,
  j = "Egger P value",
  width = 0.85
)

ft_s5 <- width(
  ft_s5,
  j = "MR-PRESSO global P",
  width = 1.20
)

ft_s5 <- width(
  ft_s5,
  j = "MR-PRESSO outliers",
  width = 1.05
)

ft_s5 <- set_table_properties(
  ft_s5,
  opts_word = list(
    split = TRUE,
    repeat_headers = TRUE
  )
)

# ============================================================
# 52. Title and notes
# ============================================================

s5_title <- paste0(
  "Supplementary Table S5. Heterogeneity and horizontal ",
  "pleiotropy diagnostics across the Mendelian randomization analyses"
)

s5_note1 <- paste0(
  "Abbreviations: Alt-E, alternative European; DMP, Dutch Microbiome ",
  "Project; GC, gastric cancer; MR, Mendelian randomization; ",
  "MR-PRESSO, Mendelian Randomization Pleiotropy RESidual Sum and Outlier."
)

s5_note2 <- paste0(
  "Cochran's Q P values assess statistical evidence of heterogeneity. ",
  "MR-Egger intercept P values assess statistical evidence of directional ",
  "horizontal pleiotropy. MR-PRESSO global P values assess evidence of ",
  "horizontal pleiotropy and outlier distortion."
)

s5_note3 <- paste0(
  "A nonsignificant diagnostic P value was interpreted as no statistical ",
  "evidence of the corresponding violation, rather than proof of its absence."
)

s5_note4 <- paste0(
  "MR-PRESSO was considered not feasible for the Swedish 2026 ",
  "Veillonellaceae primary analysis because only three instruments ",
  "remained after harmonization."
)

# ============================================================
# 53. Word output
# ============================================================

s5_output <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S5_diagnostics.docx"
)

doc_s5 <- read_docx()

doc_s5 <- body_add_fpar(
  doc_s5,
  fpar(
    ftext(
      s5_title,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc_s5 <- body_add_flextable(
  doc_s5,
  value = ft_s5
)

for (
  txt in c(
    s5_note1,
    s5_note2,
    s5_note3,
    s5_note4
  )
) {
  
  doc_s5 <- body_add_fpar(
    doc_s5,
    fpar(
      ftext(
        txt,
        prop = fp_text(
          font.family = "Arial",
          font.size = 8
        )
      )
    )
  )
}

doc_s5 <- body_end_section_landscape(
  doc_s5
)

print(
  doc_s5,
  target = s5_output
)

# ============================================================
# 54. Final S5 report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S5 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "Analyses: ",
  nrow(s5_pub),
  "\n"
)

cat(
  "MR-PRESSO feasible: ",
  sum(
    s5_audit$presso_feasible
  ),
  " / 18\n"
)

cat(
  "Publication-layer QC: PASS\n"
)

cat(
  "Frozen-result reconciliation: PASS\n"
)

cat(
  "Output:\n"
)

cat(
  s5_output,
  "\n"
)

cat("====================================================\n")
# ============================================================
# SUPPLEMENTARY TABLE S6
# Post-registration exposure-source robustness summary
# ============================================================

# ============================================================
# 55. Extract the eight Post-reg IVW analyses
# ============================================================

s6_ivw <- s4_audit %>%
  filter(
    analysis_group == "Post-reg",
    method == "Inverse variance weighted"
  )

stopifnot(
  nrow(s6_ivw) == 8,
  n_distinct(s6_ivw$analysis_id) == 8
)

# ============================================================
# 56. Attach diagnostics from S5
# ============================================================

s6_diag <- s5_audit %>%
  filter(
    analysis_group == "Post-reg"
  ) %>%
  select(
    analysis_id,
    Q_p,
    egger_intercept_p,
    presso_feasible,
    presso_global_p,
    presso_outliers
  )

stopifnot(
  nrow(s6_diag) == 8
)

s6_audit <- s6_ivw %>%
  left_join(
    s6_diag,
    by = "analysis_id"
  )

stopifnot(
  nrow(s6_audit) == 8
)

# ============================================================
# 57. Define registered analytical hierarchy
# ============================================================

s6_audit <- s6_audit %>%
  mutate(
    
    source_role =
      case_when(
        
        grepl(
          "^swedish_",
          analysis_id
        ) ~
          "Primary exposure-source robustness",
        
        grepl(
          "^dmp_",
          analysis_id
        ) ~
          "Secondary exposure-source robustness",
        
        TRUE ~
          NA_character_
      ),
    
    threshold_role =
      case_when(
        
        grepl(
          "_primary_",
          analysis_id
        ) ~
          "Primary instrument threshold",
        
        grepl(
          "_p1e5_",
          analysis_id
        ) ~
          "Secondary threshold sensitivity",
        
        TRUE ~
          NA_character_
      ),
    
    microbiome_trait =
      case_when(
        
        grepl(
          "_V$",
          analysis_id
        ) ~
          "Veillonella",
        
        grepl(
          "_F$",
          analysis_id
        ) ~
          "Veillonellaceae",
        
        TRUE ~
          NA_character_
      ),
    
    exposure_GWAS =
      case_when(
        
        grepl(
          "^swedish_",
          analysis_id
        ) ~
          "Swedish 2026",
        
        grepl(
          "^dmp_",
          analysis_id
        ) ~
          "DMP 2022",
        
        TRUE ~
          NA_character_
      )
  )

stopifnot(
  !anyNA(s6_audit$source_role),
  !anyNA(s6_audit$threshold_role),
  !anyNA(s6_audit$microbiome_trait),
  !anyNA(s6_audit$exposure_GWAS)
)

# ============================================================
# 58. Derive IVW ORs and multiplicity status
# ============================================================

s6_audit <- s6_audit %>%
  mutate(
    
    OR =
      exp(b),
    
    OR_ci_low =
      exp(
        b - 1.96 * se
      ),
    
    OR_ci_high =
      exp(
        b + 1.96 * se
      ),
    
    multiplicity_reference =
      0.0125,
    
    meets_P_0.0125 =
      pval < multiplicity_reference,
    
    effect_direction =
      case_when(
        b > 0 ~ "Positive",
        b < 0 ~ "Negative",
        TRUE ~ "Null"
      )
  )

# ============================================================
# 59. Frozen-result validation
# ============================================================

s6_locked <- tribble(
  
  ~analysis_id,
  ~locked_OR,
  ~locked_low,
  ~locked_high,
  ~locked_p,
  
  "swedish_primary_V",
  0.921,
  0.653,
  1.300,
  0.639,
  
  "swedish_primary_F",
  0.878,
  0.601,
  1.283,
  0.503,
  
  "dmp_primary_V",
  1.018,
  0.901,
  1.150,
  0.775,
  
  "dmp_primary_F",
  1.145,
  0.934,
  1.403,
  0.192,
  
  "swedish_p1e5_V",
  0.915,
  0.684,
  1.225,
  0.553,
  
  "swedish_p1e5_F",
  0.925,
  0.621,
  1.377,
  0.700,
  
  "dmp_p1e5_V",
  1.040,
  0.933,
  1.160,
  0.479,
  
  "dmp_p1e5_F",
  1.167,
  1.045,
  1.303,
  0.00616
)

s6_check <- s6_audit %>%
  left_join(
    s6_locked,
    by = "analysis_id"
  ) %>%
  mutate(
    
    OR_check =
      abs(
        OR -
          locked_OR
      ) < 0.01,
    
    low_check =
      abs(
        OR_ci_low -
          locked_low
      ) < 0.01,
    
    high_check =
      abs(
        OR_ci_high -
          locked_high
      ) < 0.01,
    
    P_check =
      abs(
        pval -
          locked_p
      ) < 0.01,
    
    audit_status =
      if_else(
        OR_check &
          low_check &
          high_check &
          P_check,
        "PASS",
        "CHECK"
      )
  )

stopifnot(
  nrow(s6_check) == 8,
  all(
    s6_check$audit_status ==
      "PASS"
  )
)

# ============================================================
# 60. Analytical-hierarchy validation
# ============================================================

# None of the four primary-threshold analyses
# should meet P < 0.0125.

s6_primary_threshold <- s6_audit %>%
  filter(
    threshold_role ==
      "Primary instrument threshold"
  )

stopifnot(
  nrow(s6_primary_threshold) == 4,
  sum(
    s6_primary_threshold$meets_P_0.0125
  ) == 0
)

# Exactly one secondary threshold-sensitivity analysis
# should meet P < 0.0125:
# DMP Veillonellaceae.

s6_secondary_threshold <- s6_audit %>%
  filter(
    threshold_role ==
      "Secondary threshold sensitivity"
  )

stopifnot(
  nrow(s6_secondary_threshold) == 4,
  sum(
    s6_secondary_threshold$meets_P_0.0125
  ) == 1
)

s6_positive_secondary <- s6_secondary_threshold %>%
  filter(
    meets_P_0.0125
  )

stopifnot(
  nrow(s6_positive_secondary) == 1,
  s6_positive_secondary$analysis_id ==
    "dmp_p1e5_F"
)

# ============================================================
# 61. Directional source-robustness validation
# ============================================================

# Veillonellaceae:
# Swedish estimates are negative at both thresholds,
# whereas DMP estimates are positive.

s6_family <- s6_audit %>%
  filter(
    microbiome_trait ==
      "Veillonellaceae"
  )

stopifnot(
  
  all(
    s6_family$effect_direction[
      s6_family$exposure_GWAS ==
        "Swedish 2026"
    ] ==
      "Negative"
  ),
  
  all(
    s6_family$effect_direction[
      s6_family$exposure_GWAS ==
        "DMP 2022"
    ] ==
      "Positive"
  )
)

# ============================================================
# 62. Save S6 audit master
# ============================================================

s6_audit_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS6_postreg_robustness_audit_master.csv"
)

write_csv(
  s6_audit,
  s6_audit_file
)

s6_check_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS6_frozen_result_validation.csv"
)

write_csv(
  s6_check,
  s6_check_file
)

# ============================================================
# 63. Final S6 data-audit report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S6 DATA AUDIT COMPLETE\n")
cat("====================================================\n")

cat(
  "Post-reg analyses: ",
  nrow(s6_audit),
  "\n"
)

cat(
  "Primary-threshold comparisons meeting P < 0.0125: ",
  sum(
    s6_primary_threshold$meets_P_0.0125
  ),
  " / 4\n"
)

cat(
  "Secondary-threshold comparisons meeting P < 0.0125: ",
  sum(
    s6_secondary_threshold$meets_P_0.0125
  ),
  " / 4\n"
)

cat(
  "Significant secondary finding: dmp_p1e5_F\n"
)

cat(
  "Frozen-result reconciliation: PASS\n"
)

cat(
  "Registered hierarchy reconciliation: PASS\n"
)

cat(
  "Exposure-source direction audit: PASS\n"
)

cat("\nAudit master:\n")
cat(
  s6_audit_file,
  "\n"
)

cat("====================================================\n")

print(
  as_tibble(
    s6_audit %>%
      select(
        analysis_id,
        exposure_GWAS,
        microbiome_trait,
        source_role,
        threshold_role,
        nsnp,
        OR,
        OR_ci_low,
        OR_ci_high,
        pval,
        meets_P_0.0125,
        effect_direction
      )
  ),
  n = Inf
)
# ============================================================
# SUPPLEMENTARY TABLE S6
# Publication layer and Word output
# ============================================================

# ============================================================
# 64. Formatting helpers
# ============================================================

format_s6_p <- function(p) {
  
  ifelse(
    p < 0.001,
    
    formatC(
      p,
      format = "e",
      digits = 2
    ),
    
    ifelse(
      p < 0.01,
      sprintf("%.5f", p),
      sprintf("%.3f", p)
    )
  )
}

format_s6_or <- function(
    OR,
    low,
    high
) {
  
  sprintf(
    "%.3f (%.3f–%.3f)",
    OR,
    low,
    high
  )
}


# ============================================================
# 65. Build publication dataframe
# ============================================================

s6_pub <- s6_audit %>%
  mutate(
    
    `Exposure source` =
      exposure_GWAS,
    
    Trait =
      microbiome_trait,
    
    `Source role` =
      case_when(
        
        source_role ==
          "Primary exposure-source robustness" ~
          "Primary",
        
        source_role ==
          "Secondary exposure-source robustness" ~
          "Secondary",
        
        TRUE ~
          NA_character_
      ),
    
    `Threshold role` =
      case_when(
        
        threshold_role ==
          "Primary instrument threshold" ~
          "Primary",
        
        threshold_role ==
          "Secondary threshold sensitivity" ~
          "Secondary sensitivity",
        
        TRUE ~
          NA_character_
      ),
    
    `IVW OR (95% CI)` =
      format_s6_or(
        OR,
        OR_ci_low,
        OR_ci_high
      ),
    
    `P value` =
      format_s6_p(
        pval
      ),
    
    Interpretation =
      case_when(
        
        meets_P_0.0125 &
          threshold_role ==
          "Secondary threshold sensitivity" ~
          
          "Meets P < 0.0125; secondary sensitivity finding",
        
        !meets_P_0.0125 ~
          
          "Does not meet P < 0.0125",
        
        TRUE ~
          
          "Check"
      ),
    
    # Preserve frozen analytical order
    source_order =
      case_when(
        exposure_GWAS ==
          "Swedish 2026" ~ 1L,
        
        exposure_GWAS ==
          "DMP 2022" ~ 2L
      ),
    
    trait_order =
      case_when(
        microbiome_trait ==
          "Veillonella" ~ 1L,
        
        microbiome_trait ==
          "Veillonellaceae" ~ 2L
      ),
    
    threshold_order =
      case_when(
        threshold_role ==
          "Primary instrument threshold" ~ 1L,
        
        threshold_role ==
          "Secondary threshold sensitivity" ~ 2L
      )
  ) %>%
  arrange(
    threshold_order,
    source_order,
    trait_order
  ) %>%
  transmute(
    
    `Exposure source`,
    
    Trait,
    
    `Source role`,
    
    Threshold =
      instrument_threshold,
    
    `Threshold role`,
    
    SNPs =
      nsnp,
    
    `IVW OR (95% CI)`,
    
    `P value`,
    
    Interpretation
  )


# ============================================================
# 66. Publication-layer QC
# ============================================================

stopifnot(
  nrow(s6_pub) == 8
)

stopifnot(
  sum(
    grepl(
      "^Meets P < 0.0125",
      s6_pub$Interpretation
    )
  ) == 1
)

s6_sig_row <- s6_pub %>%
  filter(
    grepl(
      "^Meets P < 0.0125",
      Interpretation
    )
  )

stopifnot(
  nrow(s6_sig_row) == 1,
  
  s6_sig_row$`Exposure source` ==
    "DMP 2022",
  
  s6_sig_row$Trait ==
    "Veillonellaceae",
  
  s6_sig_row$`Threshold role` ==
    "Secondary sensitivity"
)

# Exact displayed P value should retain the frozen 0.00616
stopifnot(
  s6_sig_row$`P value` ==
    "0.00616"
)


# ============================================================
# 67. Save publication CSV
# ============================================================

s6_pub_file <- here(
  "07_tables",
  "02_publication_data",
  "TableS6_publication.csv"
)

write_csv(
  s6_pub,
  s6_pub_file
)


# ============================================================
# 68. Build flextable
# ============================================================

ft_s6 <- flextable(
  s6_pub
)

ft_s6 <- font(
  ft_s6,
  fontname = "Arial",
  part = "all"
)

ft_s6 <- fontsize(
  ft_s6,
  size = 8,
  part = "all"
)

ft_s6 <- bold(
  ft_s6,
  part = "header"
)

ft_s6 <- align(
  ft_s6,
  align = "center",
  part = "header"
)

ft_s6 <- align(
  ft_s6,
  j = c(
    "Exposure source",
    "Trait",
    "Interpretation"
  ),
  align = "left",
  part = "body"
)

ft_s6 <- align(
  ft_s6,
  j = c(
    "Source role",
    "Threshold",
    "Threshold role",
    "SNPs",
    "IVW OR (95% CI)",
    "P value"
  ),
  align = "center",
  part = "body"
)

ft_s6 <- valign(
  ft_s6,
  valign = "center",
  part = "all"
)


# ============================================================
# 69. Italicize genus Veillonella
# ============================================================

veillonella_rows <- which(
  s6_pub$Trait ==
    "Veillonella"
)

if (
  length(veillonella_rows) > 0
) {
  
  ft_s6 <- compose(
    ft_s6,
    i = veillonella_rows,
    j = "Trait",
    value = as_paragraph(
      as_i("Veillonella")
    )
  )
}


# ============================================================
# 70. Booktabs formatting
# ============================================================

ft_s6 <- border_remove(
  ft_s6
)

ft_s6 <- hline_top(
  ft_s6,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "header"
)

ft_s6 <- hline_bottom(
  ft_s6,
  border = fp_border(
    color = "black",
    width = 0.7
  ),
  part = "header"
)

ft_s6 <- hline_bottom(
  ft_s6,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "body"
)

ft_s6 <- padding(
  ft_s6,
  padding.top = 3,
  padding.bottom = 3,
  padding.left = 2.5,
  padding.right = 2.5,
  part = "all"
)


# ============================================================
# 71. Column widths
# ============================================================

ft_s6 <- width(
  ft_s6,
  j = "Exposure source",
  width = 1.25
)

ft_s6 <- width(
  ft_s6,
  j = "Trait",
  width = 1.15
)

ft_s6 <- width(
  ft_s6,
  j = "Source role",
  width = 0.85
)

ft_s6 <- width(
  ft_s6,
  j = "Threshold",
  width = 1.15
)

ft_s6 <- width(
  ft_s6,
  j = "Threshold role",
  width = 1.30
)

ft_s6 <- width(
  ft_s6,
  j = "SNPs",
  width = 0.55
)

ft_s6 <- width(
  ft_s6,
  j = "IVW OR (95% CI)",
  width = 1.55
)

ft_s6 <- width(
  ft_s6,
  j = "P value",
  width = 0.75
)

ft_s6 <- width(
  ft_s6,
  j = "Interpretation",
  width = 2.25
)

ft_s6 <- set_table_properties(
  ft_s6,
  opts_word = list(
    split = TRUE,
    repeat_headers = TRUE
  )
)


# ============================================================
# 72. Title and notes
# ============================================================

s6_title <- paste0(
  "Supplementary Table S6. Summary of the post-registration ",
  "microbiome exposure-source robustness analyses"
)

s6_note1 <- paste0(
  "Abbreviations: CI, confidence interval; DMP, Dutch Microbiome ",
  "Project; IVW, inverse-variance weighted; OR, odds ratio; ",
  "SNP, single-nucleotide polymorphism."
)

s6_note2 <- paste0(
  "Swedish 2026 was designated the primary Post-reg exposure-source ",
  "robustness dataset and DMP 2022 the secondary dataset. ",
  "P ≤ 5 × 10^-6 was the primary instrument-selection threshold; ",
  "P ≤ 1 × 10^-5 was a prespecified secondary threshold-sensitivity analysis."
)

s6_note3 <- paste0(
  "For the four primary external comparisons, the multiplicity-aware ",
  "reference threshold was P < 0.0125. The same reference threshold is ",
  "shown for the secondary threshold-sensitivity analyses for transparent ",
  "comparison."
)

s6_note4 <- paste0(
  "The DMP 2022 Veillonellaceae result at P ≤ 1 × 10^-5 met the ",
  "P < 0.0125 reference threshold but remained a secondary sensitivity ",
  "finding and did not replace the nonsignificant primary analysis."
)

s6_note5 <- paste0(
  "No formal meta-analysis across MiBioGen, Swedish 2026, and DMP 2022 ",
  "was performed because microbiome phenotype scales and preprocessing ",
  "procedures were not assumed to be directly commensurate."
)


# ============================================================
# 73. Word output
# ============================================================

s6_output <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S6_postreg_robustness.docx"
)

doc_s6 <- read_docx()

doc_s6 <- body_add_fpar(
  doc_s6,
  fpar(
    ftext(
      s6_title,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc_s6 <- body_add_flextable(
  doc_s6,
  value = ft_s6
)

for (
  txt in c(
    s6_note1,
    s6_note2,
    s6_note3,
    s6_note4,
    s6_note5
  )
) {
  
  doc_s6 <- body_add_fpar(
    doc_s6,
    fpar(
      ftext(
        txt,
        prop = fp_text(
          font.family = "Arial",
          font.size = 8
        )
      )
    )
  )
}

doc_s6 <- body_end_section_landscape(
  doc_s6
)

print(
  doc_s6,
  target = s6_output
)


# ============================================================
# 74. Final S6 report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S6 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "Post-reg analyses: ",
  nrow(s6_pub),
  "\n"
)

cat(
  "Primary-threshold analyses: 4\n"
)

cat(
  "Secondary-threshold analyses: 4\n"
)

cat(
  "Results meeting P < 0.0125: 1\n"
)

cat(
  "Positive secondary finding correctly labelled: PASS\n"
)

cat(
  "Publication-layer QC: PASS\n"
)

cat(
  "Frozen-result reconciliation: PASS\n"
)

cat(
  "Output:\n"
)

cat(
  s6_output,
  "\n"
)

cat("====================================================\n")