# ============================================================
# 11_core_audit_checks.R
#
# Read-only pre-manuscript audit
# Does NOT modify frozen core analyses
# ============================================================


dir.create(
  "00_admin/audit_outputs",
  recursive = TRUE,
  showWarnings = FALSE
)

intermediate_files <- list.files(
  "02_intermediate_data",
  recursive = TRUE,
  full.names = TRUE
)

audit_candidate_files <- intermediate_files[
  grepl(
    "Veillon|IV|harmon|Finn|East|Asian|reverse|GCST90018849",
    intermediate_files,
    ignore.case = TRUE
  )
]

audit_candidate_files
# ============================================================
# Return to project root
# ============================================================

#如果前面代码包找不到文件，就用下面代码get工作路径
#getwd()

#if (basename(getwd()) == "04_R_scripts") {
  setwd("..")
#}

#getwd()

#dir.exists("02_intermediate_data")
#dir.exists("05_results")
#dir.exists("00_admin")
# ============================================================
# 5. Inspect key instrument / harmonisation objects
# ============================================================

cat("\n========================================\n")
cat("KEY AUDIT OBJECTS\n")
cat("========================================\n")


# ------------------------------------------------------------
# Helper: inspect CSV
# ------------------------------------------------------------

inspect_csv <- function(path) {
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("\n----------------------------------------\n")
  cat("FILE:", path, "\n")
  cat("ROWS:", nrow(x), "\n")
  cat("COLUMNS:\n")
  print(names(x))
  
  invisible(x)
}


# ------------------------------------------------------------
# Helper: inspect RDS
# ------------------------------------------------------------

inspect_rds <- function(path) {
  
  x <- readRDS(path)
  
  cat("\n----------------------------------------\n")
  cat("FILE:", path, "\n")
  cat("CLASS:", class(x), "\n")
  
  if (is.data.frame(x)) {
    cat("ROWS:", nrow(x), "\n")
    cat("COLUMNS:\n")
    print(names(x))
  } else {
    cat("LENGTH:", length(x), "\n")
  }
  
  invisible(x)
}


# ============================================================
# 5A. PRIMARY EXPOSURE IVs
# ============================================================

primary_v_iv <- inspect_rds(
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

primary_va_iv <- inspect_rds(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_primary_p5e6.rds"
)


# ============================================================
# 5B. PRIMARY OUTCOME EXTRACTION
# ============================================================

primary_outcome <- inspect_rds(
  "02_intermediate_data/02_outcome/GCST90018849_primary_outcome_IVs.rds"
)


# ============================================================
# 5C. PRIMARY HARMONISATION
# ============================================================

primary_v_audit <- inspect_csv(
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_harmonisation_audit.csv"
)

primary_v_ready <- inspect_rds(
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_MR_ready.rds"
)

primary_va_audit <- inspect_csv(
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_harmonisation_audit.csv"
)

primary_va_ready <- inspect_rds(
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_MR_ready.rds"
)


# ============================================================
# 5D. THRESHOLD SENSITIVITY
# ============================================================

threshold_v_iv <- inspect_rds(
  "02_intermediate_data/01_exposure/Veillonella_IV_sensitivity_p1e5.rds"
)

threshold_va_iv <- inspect_rds(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_sensitivity_p1e5.rds"
)

threshold_v_audit <- inspect_csv(
  "02_intermediate_data/04_threshold_sensitivity/Veillonella_p1e5_GCST90018849_harmonisation_audit.csv"
)

threshold_v_ready <- inspect_rds(
  "02_intermediate_data/04_threshold_sensitivity/Veillonella_p1e5_GCST90018849_MR_ready.rds"
)

threshold_va_audit <- inspect_csv(
  "02_intermediate_data/04_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_harmonisation_audit.csv"
)

threshold_va_ready <- inspect_rds(
  "02_intermediate_data/04_threshold_sensitivity/Veillonellaceae_p1e5_GCST90018849_MR_ready.rds"
)


# ============================================================
# 5E. FINNGEN OUTCOME EXTRACTION
# ============================================================

finngen_outcome <- inspect_rds(
  "02_intermediate_data/02_outcome/FinnGen/finn-b-C3_STOMACH_extracted_IVs.rds"
)


# ============================================================
# 5F. EAST ASIAN OUTCOME EXTRACTION
# ============================================================

east_asian_outcome <- inspect_rds(
  "02_intermediate_data/02_outcome/EastAsian/GCST90018629_extracted_IVs.rds"
)

east_asian_missing <- inspect_csv(
  "02_intermediate_data/02_outcome/EastAsian/GCST90018629_missing_prespecified_IVs.csv"
)


# ============================================================
# 5G. REVERSE MR
# ============================================================

reverse_iv <- inspect_rds(
  "02_intermediate_data/03_reverse_MR/GCST90018849_reverse_IVs.rds"
)

reverse_v_missing <- inspect_csv(
  "02_intermediate_data/03_reverse_MR/Veillonella_missing_GC_IVs.csv"
)

reverse_va_missing <- inspect_csv(
  "02_intermediate_data/03_reverse_MR/Veillonellaceae_missing_GC_IVs.csv"
)

reverse_v_raw <- inspect_rds(
  "02_intermediate_data/03_reverse_MR/Veillonella_reverse_outcome_raw.rds"
)

reverse_va_raw <- inspect_rds(
  "02_intermediate_data/03_reverse_MR/Veillonellaceae_reverse_outcome_raw.rds"
)
# ============================================================
# 6. Inspect harmonisation decisions
# ============================================================

cat("\n========================================\n")
cat("PRIMARY VEILLONELLA HARMONISATION\n")
cat("========================================\n")
print(primary_v_audit)

cat("\n========================================\n")
cat("PRIMARY VEILLONELLACEAE HARMONISATION\n")
cat("========================================\n")
print(primary_va_audit)

cat("\n========================================\n")
cat("THRESHOLD VEILLONELLA HARMONISATION\n")
cat("========================================\n")
print(threshold_v_audit)

cat("\n========================================\n")
cat("THRESHOLD VEILLONELLACEAE HARMONISATION\n")
cat("========================================\n")
print(threshold_va_audit)

names(primary_v_audit)
names(primary_va_audit)

#mr_keep == FALSE
# ============================================================
# 7. Build compact instrument-flow audit
#    Primary + threshold sensitivity
# ============================================================


# ------------------------------------------------------------
# Helper: obtain minimum F statistic
#
# If an F-statistic column is already present, use it.
# Otherwise calculate approximate per-SNP F = beta^2 / se^2.
# ------------------------------------------------------------

get_min_F <- function(iv) {
  
  possible_F_names <- c(
    "F",
    "F_stat",
    "Fstat",
    "F.stat",
    "f_stat",
    "F_statistics"
  )
  
  existing_F <- possible_F_names[
    possible_F_names %in% names(iv)
  ]
  
  if (length(existing_F) > 0) {
    
    x <- as.numeric(iv[[existing_F[1]]])
    
    return(
      min(x, na.rm = TRUE)
    )
  }
  
  
  if (
    all(
      c(
        "beta.exposure",
        "se.exposure"
      ) %in% names(iv)
    )
  ) {
    
    F_value <-
      (iv$beta.exposure / iv$se.exposure)^2
    
    return(
      min(F_value, na.rm = TRUE)
    )
  }
  
  
  return(NA_real_)
}


# ------------------------------------------------------------
# Helper: summarize one MR instrument flow
# ------------------------------------------------------------

summarise_instrument_flow <- function(
    analysis,
    iv,
    audit,
    ready
) {
  
  initial_n <- nrow(iv)
  
  outcome_available_n <- nrow(audit)
  
  final_n <- nrow(ready)
  
  missing_outcome_n <-
    initial_n - outcome_available_n
  
  
  excluded_harmonisation_n <-
    sum(
      audit$mr_keep == FALSE,
      na.rm = TRUE
    )
  
  
  # Palindromic and/or ambiguous variants among excluded SNPs
  pal_or_ambiguous_removed <-
    sum(
      audit$mr_keep == FALSE &
        (
          audit$palindromic == TRUE |
            audit$ambiguous == TRUE
        ),
      na.rm = TRUE
    )
  
  
  # Allele-incompatibility flag among excluded SNPs
  incompatible_removed <-
    sum(
      audit$mr_keep == FALSE &
        audit$remove == TRUE,
      na.rm = TRUE
    )
  
  
  data.frame(
    
    analysis = analysis,
    
    initial_IVs = initial_n,
    
    outcome_available_IVs = outcome_available_n,
    
    missing_in_outcome = missing_outcome_n,
    
    excluded_during_harmonisation =
      excluded_harmonisation_n,
    
    palindromic_or_ambiguous_removed =
      pal_or_ambiguous_removed,
    
    allele_incompatible_removed =
      incompatible_removed,
    
    final_MR_IVs = final_n,
    
    min_F = get_min_F(iv),
    
    stringsAsFactors = FALSE
  )
}


# ============================================================
# Primary analyses
# ============================================================

flow_primary_v <- summarise_instrument_flow(
  analysis = "Primary: Veillonella -> GCST90018849",
  iv = primary_v_iv,
  audit = primary_v_audit,
  ready = primary_v_ready
)

flow_primary_va <- summarise_instrument_flow(
  analysis = "Primary: Veillonellaceae -> GCST90018849",
  iv = primary_va_iv,
  audit = primary_va_audit,
  ready = primary_va_ready
)


# ============================================================
# Threshold sensitivity analyses
# ============================================================

flow_threshold_v <- summarise_instrument_flow(
  analysis = "Threshold P<=1e-5: Veillonella -> GCST90018849",
  iv = threshold_v_iv,
  audit = threshold_v_audit,
  ready = threshold_v_ready
)

flow_threshold_va <- summarise_instrument_flow(
  analysis = "Threshold P<=1e-5: Veillonellaceae -> GCST90018849",
  iv = threshold_va_iv,
  audit = threshold_va_audit,
  ready = threshold_va_ready
)


# ============================================================
# Combine
# ============================================================

instrument_flow_forward_core <- rbind(
  flow_primary_v,
  flow_primary_va,
  flow_threshold_v,
  flow_threshold_va
)


cat("\n========================================\n")
cat("PRIMARY + THRESHOLD INSTRUMENT FLOW\n")
cat("========================================\n")

print(
  instrument_flow_forward_core,
  row.names = FALSE
)


# ============================================================
# Save audit output
# ============================================================

write.csv(
  instrument_flow_forward_core,
  "00_admin/audit_outputs/instrument_flow_primary_threshold.csv",
  row.names = FALSE
)
# ============================================================
# Fix audit output directory
# ============================================================

getwd()

dir.create(
  "00_admin/audit_outputs",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.exists("00_admin/audit_outputs")
#重新保存
write.csv(
  instrument_flow_forward_core,
  "00_admin/audit_outputs/instrument_flow_primary_threshold.csv",
  row.names = FALSE
)
#检查
file.exists(
  "00_admin/audit_outputs/instrument_flow_primary_threshold.csv"
)
#压缩表
compact_forward_flow <- instrument_flow_forward_core[
  ,
  c(
    "analysis",
    "initial_IVs",
    "outcome_available_IVs",
    "missing_in_outcome",
    "excluded_during_harmonisation",
    "palindromic_or_ambiguous_removed",
    "final_MR_IVs",
    "min_F"
  )
]

print(
  compact_forward_flow,
  row.names = FALSE
)
#检查数学闭环
with(
  instrument_flow_forward_core,
  outcome_available_IVs -
    excluded_during_harmonisation ==
    final_MR_IVs
)
# ============================================================
# 8. Locate FinnGen, East Asian and reverse audit/result files
# ============================================================

remaining_audit_files <- list.files(
  c(
    "02_intermediate_data",
    "05_results"
  ),
  recursive = TRUE,
  full.names = TRUE
)

remaining_audit_files <- remaining_audit_files[
  grepl(
    "Finn|finn|East|Asian|18629|reverse|Reverse",
    remaining_audit_files
  )
]

cat("\n========================================\n")
cat("FINNGEN / EAST ASIAN / REVERSE FILES\n")
cat("========================================\n")

remaining_audit_files
#筛选
cat("\n===== FINNGEN =====\n")
remaining_audit_files[
  grepl(
    "Finn|finn",
    remaining_audit_files
  )
]

cat("\n===== EAST ASIAN =====\n")
remaining_audit_files[
  grepl(
    "East|Asian|18629",
    remaining_audit_files
  )
]

cat("\n===== REVERSE =====\n")
remaining_audit_files[
  grepl(
    "reverse|Reverse",
    remaining_audit_files
  )
]
# ============================================================
# 9. Compact file audit for FinnGen / East Asian / Reverse
# ============================================================

cat("\n===== FINNGEN =====\n")

finngen_files <- remaining_audit_files[
  grepl(
    "Finn|finn",
    remaining_audit_files
  )
]

cat(
  paste(finngen_files, collapse = "\n"),
  "\n"
)


cat("\n===== EAST ASIAN =====\n")

east_asian_files <- remaining_audit_files[
  grepl(
    "East|Asian|18629",
    remaining_audit_files
  )
]

cat(
  paste(east_asian_files, collapse = "\n"),
  "\n"
)


cat("\n===== REVERSE =====\n")

reverse_files <- remaining_audit_files[
  grepl(
    "reverse|Reverse",
    remaining_audit_files
  )
]

cat(
  paste(reverse_files, collapse = "\n"),
  "\n"
)
# ============================================================
# 10. Complete instrument-flow audit:
#     FinnGen + East Asian + Reverse MR
# ============================================================


# ------------------------------------------------------------
# Helper: count unique SNPs safely
# ------------------------------------------------------------

count_snps <- function(x) {
  
  if (!is.data.frame(x)) {
    stop("Object is not a data.frame.")
  }
  
  snp_col <- names(x)[
    tolower(names(x)) == "snp"
  ]
  
  if (length(snp_col) == 0) {
    stop("No SNP column found.")
  }
  
  length(
    unique(
      x[[snp_col[1]]]
    )
  )
}


# ------------------------------------------------------------
# Helper: get SNP vector
# ------------------------------------------------------------

get_snps <- function(x) {
  
  snp_col <- names(x)[
    tolower(names(x)) == "snp"
  ]
  
  if (length(snp_col) == 0) {
    stop("No SNP column found.")
  }
  
  unique(
    x[[snp_col[1]]]
  )
}


# ------------------------------------------------------------
# Helper: locate exactly one MR result file
# Avoid substring collision:
# Veillonella must NOT match Veillonellaceae
# ------------------------------------------------------------

find_mr_file <- function(
    directory,
    taxon
) {
  
  x <- list.files(
    directory,
    pattern = "_MR\\.csv$",
    full.names = TRUE
  )
  
  # Match taxon as a complete filename token
  pattern <- paste0(
    "(^|_)",
    taxon,
    "(_|\\.)"
  )
  
  x <- x[
    grepl(
      pattern,
      basename(x),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]
  
  if (length(x) != 1) {
    
    cat("\nFiles found for:", taxon, "\n")
    print(x)
    
    stop(
      "Expected exactly one MR result file."
    )
  }
  
  x
}
#测试
find_mr_file(
  "05_results/03_finngen_sensitivity",
  "Veillonella"
)

find_mr_file(
  "05_results/03_finngen_sensitivity",
  "Veillonellaceae"
)
# ------------------------------------------------------------
# Helper: obtain IVW nsnp from saved MR result
# ------------------------------------------------------------

get_ivw_nsnp <- function(path) {
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  ivw <- x[
    x$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]
  
  if (nrow(ivw) != 1) {
    stop(
      paste(
        "Could not identify exactly one IVW row in:",
        path
      )
    )
  }
  
  as.integer(ivw$nsnp)
}


# ============================================================
# 10A. FINNGEN
# ============================================================

cat("\n========================================\n")
cat("FINNGEN INSTRUMENT FLOW\n")
cat("========================================\n")


# Exposure IV sets
v_primary_snps <-
  get_snps(primary_v_iv)

va_primary_snps <-
  get_snps(primary_va_iv)


# SNPs actually available in FinnGen outcome
finngen_snps <-
  get_snps(finngen_outcome)


finngen_v_available <-
  length(
    intersect(
      v_primary_snps,
      finngen_snps
    )
  )

finngen_va_available <-
  length(
    intersect(
      va_primary_snps,
      finngen_snps
    )
  )


# Locate final MR result files
finngen_v_mr_file <- find_mr_file(
  "05_results/03_finngen_sensitivity",
  "Veillonella"
)

finngen_va_mr_file <- find_mr_file(
  "05_results/03_finngen_sensitivity",
  "Veillonellaceae"
)


# Final IV counts from IVW result
finngen_v_final <-
  get_ivw_nsnp(finngen_v_mr_file)

finngen_va_final <-
  get_ivw_nsnp(finngen_va_mr_file)


flow_finngen_v <- data.frame(
  
  analysis =
    "FinnGen sensitivity: Veillonella -> gastric cancer",
  
  initial_IVs =
    nrow(primary_v_iv),
  
  outcome_available_IVs =
    finngen_v_available,
  
  missing_in_outcome =
    nrow(primary_v_iv) -
    finngen_v_available,
  
  excluded_during_harmonisation =
    finngen_v_available -
    finngen_v_final,
  
  final_MR_IVs =
    finngen_v_final,
  
  min_F =
    get_min_F(primary_v_iv),
  
  stringsAsFactors = FALSE
)


flow_finngen_va <- data.frame(
  
  analysis =
    "FinnGen sensitivity: Veillonellaceae -> gastric cancer",
  
  initial_IVs =
    nrow(primary_va_iv),
  
  outcome_available_IVs =
    finngen_va_available,
  
  missing_in_outcome =
    nrow(primary_va_iv) -
    finngen_va_available,
  
  excluded_during_harmonisation =
    finngen_va_available -
    finngen_va_final,
  
  final_MR_IVs =
    finngen_va_final,
  
  min_F =
    get_min_F(primary_va_iv),
  
  stringsAsFactors = FALSE
)



# ============================================================
# 10B. EAST ASIAN
# ============================================================

cat("\n========================================\n")
cat("EAST ASIAN INSTRUMENT FLOW\n")
cat("========================================\n")


east_asian_snps <-
  get_snps(east_asian_outcome)


east_v_available <-
  length(
    intersect(
      v_primary_snps,
      east_asian_snps
    )
  )

east_va_available <-
  length(
    intersect(
      va_primary_snps,
      east_asian_snps
    )
  )


east_v_mr_file <- find_mr_file(
  "05_results/04_east_asian_sensitivity",
  "Veillonella"
)

east_va_mr_file <- find_mr_file(
  "05_results/04_east_asian_sensitivity",
  "Veillonellaceae"
)


east_v_final <-
  get_ivw_nsnp(east_v_mr_file)

east_va_final <-
  get_ivw_nsnp(east_va_mr_file)


flow_east_v <- data.frame(
  
  analysis =
    "East Asian sensitivity: Veillonella -> gastric cancer",
  
  initial_IVs =
    nrow(primary_v_iv),
  
  outcome_available_IVs =
    east_v_available,
  
  missing_in_outcome =
    nrow(primary_v_iv) -
    east_v_available,
  
  excluded_during_harmonisation =
    east_v_available -
    east_v_final,
  
  final_MR_IVs =
    east_v_final,
  
  min_F =
    get_min_F(primary_v_iv),
  
  stringsAsFactors = FALSE
)


flow_east_va <- data.frame(
  
  analysis =
    "East Asian sensitivity: Veillonellaceae -> gastric cancer",
  
  initial_IVs =
    nrow(primary_va_iv),
  
  outcome_available_IVs =
    east_va_available,
  
  missing_in_outcome =
    nrow(primary_va_iv) -
    east_va_available,
  
  excluded_during_harmonisation =
    east_va_available -
    east_va_final,
  
  final_MR_IVs =
    east_va_final,
  
  min_F =
    get_min_F(primary_va_iv),
  
  stringsAsFactors = FALSE
)



# ============================================================
# 10C. REVERSE MR
# ============================================================

cat("\n========================================\n")
cat("REVERSE MR INSTRUMENT FLOW\n")
cat("========================================\n")


reverse_initial_n <-
  nrow(reverse_iv)


reverse_v_available <-
  nrow(reverse_v_raw)

reverse_va_available <-
  nrow(reverse_va_raw)


reverse_v_missing_n <-
  nrow(reverse_v_missing)

reverse_va_missing_n <-
  nrow(reverse_va_missing)


reverse_v_mr_file <- find_mr_file(
  "05_results/05_reverse_MR",
  "GC_to_Veillonella"
)

reverse_va_mr_file <- find_mr_file(
  "05_results/05_reverse_MR",
  "GC_to_Veillonellaceae"
)


reverse_v_final <-
  get_ivw_nsnp(reverse_v_mr_file)

reverse_va_final <-
  get_ivw_nsnp(reverse_va_mr_file)


flow_reverse_v <- data.frame(
  
  analysis =
    "Reverse: gastric cancer -> Veillonella",
  
  initial_IVs =
    reverse_initial_n,
  
  outcome_available_IVs =
    reverse_v_available,
  
  missing_in_outcome =
    reverse_v_missing_n,
  
  excluded_during_harmonisation =
    reverse_v_available -
    reverse_v_final,
  
  final_MR_IVs =
    reverse_v_final,
  
  min_F =
    get_min_F(reverse_iv),
  
  stringsAsFactors = FALSE
)


flow_reverse_va <- data.frame(
  
  analysis =
    "Reverse: gastric cancer -> Veillonellaceae",
  
  initial_IVs =
    reverse_initial_n,
  
  outcome_available_IVs =
    reverse_va_available,
  
  missing_in_outcome =
    reverse_va_missing_n,
  
  excluded_during_harmonisation =
    reverse_va_available -
    reverse_va_final,
  
  final_MR_IVs =
    reverse_va_final,
  
  min_F =
    get_min_F(reverse_iv),
  
  stringsAsFactors = FALSE
)



# ============================================================
# 11. Combine COMPLETE instrument-flow audit
# ============================================================

# First simplify the already completed forward table
forward_core_simple <-
  instrument_flow_forward_core[
    ,
    c(
      "analysis",
      "initial_IVs",
      "outcome_available_IVs",
      "missing_in_outcome",
      "excluded_during_harmonisation",
      "final_MR_IVs",
      "min_F"
    )
  ]


instrument_flow_complete <- rbind(
  
  forward_core_simple,
  
  flow_finngen_v,
  flow_finngen_va,
  
  flow_east_v,
  flow_east_va,
  
  flow_reverse_v,
  flow_reverse_va
)


# ============================================================
# 12. Mathematical consistency check
# ============================================================

instrument_flow_complete$flow_consistent <-
  with(
    instrument_flow_complete,
    
    initial_IVs -
      missing_in_outcome -
      excluded_during_harmonisation ==
      final_MR_IVs
  )


cat("\n========================================\n")
cat("COMPLETE CORE INSTRUMENT FLOW\n")
cat("========================================\n")

print(
  instrument_flow_complete,
  row.names = FALSE
)


cat("\n========================================\n")
cat("FLOW CONSISTENCY\n")
cat("========================================\n")

print(
  instrument_flow_complete[
    ,
    c(
      "analysis",
      "initial_IVs",
      "outcome_available_IVs",
      "missing_in_outcome",
      "excluded_during_harmonisation",
      "final_MR_IVs",
      "min_F",
      "flow_consistent"
    )
  ],
  row.names = FALSE
)



# ============================================================
# 13. Save complete audit table
# ============================================================

dir.create(
  "00_admin/audit_outputs",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  instrument_flow_complete,
  "00_admin/audit_outputs/instrument_flow_audit_complete.csv",
  row.names = FALSE
)


cat(
  "\nSaved:\n",
  "00_admin/audit_outputs/instrument_flow_audit_complete.csv\n"
)


file.exists(
  "00_admin/audit_outputs/instrument_flow_audit_complete.csv"
)

#检查
# ============================================================
# Final checks for Instrument / Harmonisation Audit
# ============================================================

nrow(instrument_flow_complete)

all(
  instrument_flow_complete$flow_consistent
)

all(
  instrument_flow_complete$min_F > 10
)

instrument_flow_complete[
  ,
  c(
    "analysis",
    "initial_IVs",
    "outcome_available_IVs",
    "missing_in_outcome",
    "excluded_during_harmonisation",
    "final_MR_IVs",
    "min_F",
    "flow_consistent"
  )
]
# ============================================================
# 14. RESULT PROVENANCE AUDIT
#     Raw MR result files -> all_MR_results -> IVW_summary
# ============================================================


# ============================================================
# 14A. Locate the 10 source MR result files
# ============================================================

mr_source_dirs <- c(
  "05_results/01_primary_MR",
  "05_results/02_threshold_sensitivity",
  "05_results/03_finngen_sensitivity",
  "05_results/04_east_asian_sensitivity",
  "05_results/05_reverse_MR"
)

mr_source_files <- unlist(
  lapply(
    mr_source_dirs,
    function(d) {
      list.files(
        d,
        pattern = "_MR\\.csv$",
        full.names = TRUE
      )
    }
  )
)

cat("\n========================================\n")
cat("SOURCE MR FILES\n")
cat("========================================\n")

print(mr_source_files)

cat("\nNumber of source MR files:\n")
print(length(mr_source_files))
# ============================================================
# 14B. Extract IVW row from every source MR file
# ============================================================

extract_source_ivw <- function(path) {
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  ivw <- x[
    x$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]
  
  if (nrow(ivw) != 1) {
    
    stop(
      paste(
        "Expected exactly one IVW row in:",
        path
      )
    )
  }
  
  
  # Optional OR columns
  get_optional <- function(name) {
    
    if (name %in% names(ivw)) {
      return(as.numeric(ivw[[name]][1]))
    }
    
    NA_real_
  }
  
  
  data.frame(
    
    source_file =
      basename(path),
    
    source_folder =
      basename(dirname(path)),
    
    nsnp =
      as.integer(ivw$nsnp[1]),
    
    b =
      as.numeric(ivw$b[1]),
    
    se =
      as.numeric(ivw$se[1]),
    
    pval =
      as.numeric(ivw$pval[1]),
    
    or =
      get_optional("or"),
    
    or_lci95 =
      get_optional("or_lci95"),
    
    or_uci95 =
      get_optional("or_uci95"),
    
    stringsAsFactors = FALSE
  )
}


source_ivw <- do.call(
  rbind,
  lapply(
    mr_source_files,
    extract_source_ivw
  )
)


cat("\n========================================\n")
cat("IVW VALUES FROM SOURCE FILES\n")
cat("========================================\n")

print(
  source_ivw[
    ,
    c(
      "source_folder",
      "source_file",
      "nsnp",
      "b",
      "se",
      "pval"
    )
  ],
  row.names = FALSE
)
# ============================================================
# 14C. Read frozen summary files
# ============================================================

all_mr_saved <- read.csv(
  "05_results/06_summary/all_MR_results.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ivw_summary_saved <- read.csv(
  "05_results/06_summary/IVW_summary.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


all_mr_ivw_saved <- all_mr_saved[
  all_mr_saved$method == "Inverse variance weighted",
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("SUMMARY ROW COUNTS\n")
cat("========================================\n")

cat(
  "Source MR files:",
  nrow(source_ivw),
  "\n"
)

cat(
  "IVW rows in all_MR_results:",
  nrow(all_mr_ivw_saved),
  "\n"
)

cat(
  "Rows in IVW_summary:",
  nrow(ivw_summary_saved),
  "\n"
)
# ============================================================
# 14D. Match every raw IVW result to all_MR_results
# ============================================================

tol <- 1e-10


match_one_source_result <- function(i) {
  
  s <- source_ivw[i, ]
  
  
  candidates <-
    which(
      
      all_mr_ivw_saved$nsnp == s$nsnp &
        
        abs(
          all_mr_ivw_saved$b -
            s$b
        ) < tol &
        
        abs(
          all_mr_ivw_saved$se -
            s$se
        ) < tol &
        
        abs(
          all_mr_ivw_saved$pval -
            s$pval
        ) < tol
    )
  
  
  if (length(candidates) == 1) {
    
    j <- candidates[1]
    
    data.frame(
      
      source_folder =
        s$source_folder,
      
      source_file =
        s$source_file,
      
      matched_rows =
        1,
      
      matched_analysis =
        if ("analysis" %in% names(all_mr_ivw_saved))
          all_mr_ivw_saved$analysis[j]
      else
        NA,
      
      matched_exposure =
        if ("exposure_label" %in% names(all_mr_ivw_saved))
          all_mr_ivw_saved$exposure_label[j]
      else
        NA,
      
      matched_outcome =
        if ("outcome_label" %in% names(all_mr_ivw_saved))
          all_mr_ivw_saved$outcome_label[j]
      else
        NA,
      
      nsnp =
        s$nsnp,
      
      b =
        s$b,
      
      se =
        s$se,
      
      pval =
        s$pval,
      
      exact_unique_match =
        TRUE,
      
      stringsAsFactors = FALSE
    )
    
  } else {
    
    data.frame(
      
      source_folder =
        s$source_folder,
      
      source_file =
        s$source_file,
      
      matched_rows =
        length(candidates),
      
      matched_analysis =
        NA,
      
      matched_exposure =
        NA,
      
      matched_outcome =
        NA,
      
      nsnp =
        s$nsnp,
      
      b =
        s$b,
      
      se =
        s$se,
      
      pval =
        s$pval,
      
      exact_unique_match =
        FALSE,
      
      stringsAsFactors = FALSE
    )
  }
}


provenance_raw_to_allMR <- do.call(
  rbind,
  lapply(
    seq_len(nrow(source_ivw)),
    match_one_source_result
  )
)


cat("\n========================================\n")
cat("RAW FILE -> ALL_MR_RESULTS MATCHING\n")
cat("========================================\n")

print(
  provenance_raw_to_allMR[
    ,
    c(
      "source_folder",
      "source_file",
      "matched_analysis",
      "nsnp",
      "b",
      "pval",
      "exact_unique_match"
    )
  ],
  row.names = FALSE
)
# ============================================================
# 14E. Global provenance checks
# ============================================================

check_source_file_count <-
  nrow(source_ivw) == 10

check_allMR_ivw_count <-
  nrow(all_mr_ivw_saved) == 10

check_IVW_summary_count <-
  nrow(ivw_summary_saved) == 10

check_unique_raw_matches <-
  all(
    provenance_raw_to_allMR$exact_unique_match
  )


cat("\n========================================\n")
cat("RESULT PROVENANCE CHECKS\n")
cat("========================================\n")

cat(
  "10 source MR files: ",
  check_source_file_count,
  "\n"
)

cat(
  "10 IVW rows in all_MR_results: ",
  check_allMR_ivw_count,
  "\n"
)

cat(
  "10 rows in IVW_summary: ",
  check_IVW_summary_count,
  "\n"
)

cat(
  "Every source IVW uniquely matches all_MR_results: ",
  check_unique_raw_matches,
  "\n"
)
#验证
names(ivw_summary_saved)
# ============================================================
# 14F. Correct IVW_summary provenance check
#     Use unique analysis + exposure + outcome keys
# ============================================================

key_vars <- c(
  "analysis",
  "exposure_label",
  "outcome_label"
)

stopifnot(
  all(key_vars %in% names(ivw_summary_saved)),
  all(key_vars %in% names(all_mr_ivw_saved))
)


# Confirm that the key is unique in both tables
summary_key <-
  interaction(
    ivw_summary_saved[, key_vars],
    drop = TRUE
  )

raw_key <-
  interaction(
    all_mr_ivw_saved[, key_vars],
    drop = TRUE
  )


cat("\nDuplicate keys in IVW_summary:\n")
print(anyDuplicated(summary_key))

cat("\nDuplicate keys in all_MR_results IVW rows:\n")
print(anyDuplicated(raw_key))

#合并
provenance_summary_check <- merge(
  
  ivw_summary_saved,
  
  all_mr_ivw_saved,
  
  by = key_vars,
  
  suffixes = c(
    ".summary",
    ".raw"
  ),
  
  all.x = TRUE
)
#检查
nrow(provenance_summary_check)
#再次核对nsnp和p
provenance_summary_check$nsnp_match <-
  provenance_summary_check$nsnp.summary ==
  provenance_summary_check$nsnp.raw


provenance_summary_check$pval_match <-
  abs(
    provenance_summary_check$pval.summary -
      provenance_summary_check$pval.raw
  ) < 1e-10


print(
  provenance_summary_check[
    ,
    c(
      "analysis",
      "exposure_label",
      "outcome_label",
      "nsnp.summary",
      "nsnp.raw",
      "pval.summary",
      "pval.raw",
      "nsnp_match",
      "pval_match"
    )
  ],
  row.names = FALSE
)


cat(
  "\nAll nsnp consistent: ",
  all(provenance_summary_check$nsnp_match),
  "\n"
)

cat(
  "All P values consistent: ",
  all(provenance_summary_check$pval_match),
  "\n"
)
# ============================================================
# 14G. Robust effect-estimate provenance check
#      Avoid merge suffix/name ambiguity
# ============================================================

key_vars <- c(
  "analysis",
  "exposure_label",
  "outcome_label"
)

# ------------------------------------------------------------
# Create unique keys
# ------------------------------------------------------------

summary_key <- do.call(
  paste,
  c(
    ivw_summary_saved[key_vars],
    sep = "|||"
  )
)

raw_key <- do.call(
  paste,
  c(
    all_mr_ivw_saved[key_vars],
    sep = "|||"
  )
)

# Each summary row should map to exactly one raw IVW row
raw_index <- match(
  summary_key,
  raw_key
)

cat(
  "All IVW_summary rows matched to raw IVW rows: ",
  all(!is.na(raw_index)),
  "\n"
)

cat(
  "Number of matched rows: ",
  sum(!is.na(raw_index)),
  "\n"
)
# ============================================================
# Pull raw beta / SE by unique key
# ============================================================

ivw_summary_check <- ivw_summary_saved

ivw_summary_check$b_raw <-
  as.numeric(
    all_mr_ivw_saved$b[raw_index]
  )

ivw_summary_check$se_raw <-
  as.numeric(
    all_mr_ivw_saved$se[raw_index]
  )

ivw_summary_check$pval_raw <-
  as.numeric(
    all_mr_ivw_saved$pval[raw_index]
  )

ivw_summary_check$nsnp_raw <-
  as.integer(
    all_mr_ivw_saved$nsnp[raw_index]
  )

sapply(
  ivw_summary_check[
    ,
    c(
      "b_raw",
      "se_raw",
      "pval_raw",
      "nsnp_raw"
    )
  ],
  function(x) sum(is.na(x))
)
# ============================================================
# Verify reported estimate
# ============================================================

ivw_summary_check$expected_estimate <-
  ifelse(
    ivw_summary_check$direction == "Forward",
    
    exp(
      ivw_summary_check$b_raw
    ),
    
    ivw_summary_check$b_raw
  )


ivw_summary_check$estimate_match <-
  abs(
    as.numeric(ivw_summary_check$estimate) -
      ivw_summary_check$expected_estimate
  ) < 1e-8


cat(
  "All reported effect estimates consistent: ",
  all(ivw_summary_check$estimate_match),
  "\n"
)
# ============================================================
# Verify reported 95% confidence intervals
# ============================================================

ivw_summary_check$expected_ci_lower <-
  ifelse(
    ivw_summary_check$direction == "Forward",
    
    exp(
      ivw_summary_check$b_raw -
        1.96 * ivw_summary_check$se_raw
    ),
    
    ivw_summary_check$b_raw -
      1.96 * ivw_summary_check$se_raw
  )


ivw_summary_check$expected_ci_upper <-
  ifelse(
    ivw_summary_check$direction == "Forward",
    
    exp(
      ivw_summary_check$b_raw +
        1.96 * ivw_summary_check$se_raw
    ),
    
    ivw_summary_check$b_raw +
      1.96 * ivw_summary_check$se_raw
  )


ivw_summary_check$ci_lower_match <-
  abs(
    as.numeric(ivw_summary_check$ci_lower) -
      ivw_summary_check$expected_ci_lower
  ) < 1e-6


ivw_summary_check$ci_upper_match <-
  abs(
    as.numeric(ivw_summary_check$ci_upper) -
      ivw_summary_check$expected_ci_upper
  ) < 1e-6


cat(
  "All lower CIs consistent: ",
  all(ivw_summary_check$ci_lower_match),
  "\n"
)

cat(
  "All upper CIs consistent: ",
  all(ivw_summary_check$ci_upper_match),
  "\n"
)

ivw_summary_check$nsnp_match <-
  as.integer(ivw_summary_check$nsnp) ==
  ivw_summary_check$nsnp_raw


ivw_summary_check$pval_match <-
  abs(
    as.numeric(ivw_summary_check$pval) -
      ivw_summary_check$pval_raw
  ) < 1e-10


cat(
  "All nsnp consistent: ",
  all(ivw_summary_check$nsnp_match),
  "\n"
)

cat(
  "All P values consistent: ",
  all(ivw_summary_check$pval_match),
  "\n"
)
# ============================================================
# 15. QC RESULT PROVENANCE AUDIT
#
# Source heterogeneity / Egger / MR-PRESSO files
#                   ↓
#             QC_summary.csv
#
# Read-only audit
# ============================================================


# ============================================================
# 15A. Locate all source heterogeneity files
# ============================================================

qc_source_dirs <- c(
  "05_results/01_primary_MR",
  "05_results/02_threshold_sensitivity",
  "05_results/03_finngen_sensitivity",
  "05_results/04_east_asian_sensitivity",
  "05_results/05_reverse_MR"
)

heterogeneity_files <- unlist(
  lapply(
    qc_source_dirs,
    function(d) {
      list.files(
        d,
        pattern = "_heterogeneity\\.csv$",
        full.names = TRUE
      )
    }
  )
)

cat("\n========================================\n")
cat("SOURCE HETEROGENEITY FILES\n")
cat("========================================\n")

print(heterogeneity_files)

cat("\nNumber of heterogeneity files:\n")
print(length(heterogeneity_files))
# ============================================================
# 15B. Build QC source-file triplets
# ============================================================

qc_triplets <- lapply(
  heterogeneity_files,
  function(het_path) {
    
    folder <- dirname(het_path)
    
    prefix <- sub(
      "_heterogeneity\\.csv$",
      "",
      basename(het_path)
    )
    
    egger_path <- file.path(
      folder,
      paste0(
        prefix,
        "_Egger_intercept.csv"
      )
    )
    
    presso_path <- file.path(
      folder,
      paste0(
        prefix,
        "_MRPRESSO.rds"
      )
    )
    
    data.frame(
      folder = folder,
      prefix = prefix,
      heterogeneity_file = het_path,
      egger_file = egger_path,
      presso_file = presso_path,
      egger_exists = file.exists(egger_path),
      presso_exists = file.exists(presso_path),
      stringsAsFactors = FALSE
    )
  }
)

qc_triplets <- do.call(
  rbind,
  qc_triplets
)


cat("\n========================================\n")
cat("QC FILE TRIPLETS\n")
cat("========================================\n")

print(
  qc_triplets[
    ,
    c(
      "prefix",
      "egger_exists",
      "presso_exists"
    )
  ],
  row.names = FALSE
)

#检查
all(qc_triplets$egger_exists)
all(qc_triplets$presso_exists)
# ============================================================
# 15C. Extract QC values from source files
# ============================================================

extract_qc_source <- function(
    het_path,
    egger_path,
    presso_path
) {
  
  # ----------------------------------------------------------
  # Heterogeneity
  # ----------------------------------------------------------
  
  het <- read.csv(
    het_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  ivw_het <- het[
    het$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]
  
  if (nrow(ivw_het) != 1) {
    stop(
      paste(
        "Expected exactly one IVW heterogeneity row:",
        het_path
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # MR-Egger intercept
  # ----------------------------------------------------------
  
  egger <- read.csv(
    egger_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (nrow(egger) != 1) {
    stop(
      paste(
        "Expected exactly one Egger intercept row:",
        egger_path
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # MR-PRESSO
  # ----------------------------------------------------------
  
  presso <- readRDS(
    presso_path
  )
  
  global_p <-
    presso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  
  global_p <- as.numeric(
    sub(
      "^<",
      "",
      as.character(global_p)
    )
  )
  
  
  # ----------------------------------------------------------
  # Return compact row
  # ----------------------------------------------------------
  
  data.frame(
    
    source_folder =
      basename(dirname(het_path)),
    
    source_prefix =
      sub(
        "_heterogeneity\\.csv$",
        "",
        basename(het_path)
      ),
    
    IVW_Q =
      as.numeric(ivw_het$Q[1]),
    
    IVW_Q_df =
      as.numeric(ivw_het$Q_df[1]),
    
    IVW_Q_p =
      as.numeric(ivw_het$Q_pval[1]),
    
    Egger_intercept =
      as.numeric(egger$egger_intercept[1]),
    
    Egger_intercept_SE =
      as.numeric(egger$se[1]),
    
    Egger_intercept_p =
      as.numeric(egger$pval[1]),
    
    MRPRESSO_global_p =
      global_p,
    
    stringsAsFactors = FALSE
  )
}
source_qc <- do.call(
  rbind,
  lapply(
    seq_len(nrow(qc_triplets)),
    function(i) {
      
      extract_qc_source(
        qc_triplets$heterogeneity_file[i],
        qc_triplets$egger_file[i],
        qc_triplets$presso_file[i]
      )
    }
  )
)


cat("\n========================================\n")
cat("QC VALUES FROM SOURCE FILES\n")
cat("========================================\n")

print(
  source_qc,
  row.names = FALSE
)
# ============================================================
# 15D. Read saved QC summary
# ============================================================

qc_summary_saved <- read.csv(
  "05_results/06_summary/QC_summary.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("QC ROW COUNTS\n")
cat("========================================\n")

cat(
  "Source QC analyses:",
  nrow(source_qc),
  "\n"
)

cat(
  "Rows in QC_summary:",
  nrow(qc_summary_saved),
  "\n"
)
# ============================================================
# 15E. Numeric provenance matching
# ============================================================

tol_qc <- 1e-8


match_one_qc <- function(i) {
  
  s <- source_qc[i, ]
  
  
  candidates <- which(
    
    abs(
      qc_summary_saved$IVW_Q -
        s$IVW_Q
    ) < tol_qc &
      
      abs(
        qc_summary_saved$IVW_Q_df -
          s$IVW_Q_df
      ) < tol_qc &
      
      abs(
        qc_summary_saved$IVW_Q_p -
          s$IVW_Q_p
      ) < tol_qc &
      
      abs(
        qc_summary_saved$Egger_intercept -
          s$Egger_intercept
      ) < tol_qc &
      
      abs(
        qc_summary_saved$Egger_intercept_SE -
          s$Egger_intercept_SE
      ) < tol_qc &
      
      abs(
        qc_summary_saved$Egger_intercept_p -
          s$Egger_intercept_p
      ) < tol_qc &
      
      abs(
        qc_summary_saved$MRPRESSO_global_p -
          s$MRPRESSO_global_p
      ) < tol_qc
  )
  
  
  data.frame(
    
    source_folder =
      s$source_folder,
    
    source_prefix =
      s$source_prefix,
    
    matched_rows =
      length(candidates),
    
    matched_analysis =
      if (length(candidates) == 1)
        qc_summary_saved$analysis[candidates]
    else
      NA,
    
    matched_comparison =
      if (
        length(candidates) == 1 &&
        "comparison" %in% names(qc_summary_saved)
      )
        qc_summary_saved$comparison[candidates]
    else
      NA,
    
    exact_unique_match =
      length(candidates) == 1,
    
    stringsAsFactors = FALSE
  )
}


qc_provenance_match <- do.call(
  rbind,
  lapply(
    seq_len(nrow(source_qc)),
    match_one_qc
  )
)
#打印
cat("\n========================================\n")
cat("SOURCE QC -> QC_SUMMARY MATCHING\n")
cat("========================================\n")

print(
  qc_provenance_match,
  row.names = FALSE
)
# ============================================================
# 15F. Global QC provenance checks
# ============================================================

check_qc_source_count <-
  nrow(source_qc) == 10

check_qc_summary_count <-
  nrow(qc_summary_saved) == 10

check_qc_unique_matches <-
  all(
    qc_provenance_match$exact_unique_match
  )


cat("\n========================================\n")
cat("QC PROVENANCE CHECKS\n")
cat("========================================\n")

cat(
  "10 source QC analyses: ",
  check_qc_source_count,
  "\n"
)

cat(
  "10 rows in QC_summary: ",
  check_qc_summary_count,
  "\n"
)

cat(
  "Every source QC result uniquely matches QC_summary: ",
  check_qc_unique_matches,
  "\n"
)
# ============================================================
# 15G. Verify derived QC flags
# ============================================================

expected_heterogeneity <-
  qc_summary_saved$IVW_Q_p < 0.05

expected_egger_pleiotropy <-
  qc_summary_saved$Egger_intercept_p < 0.05

expected_presso <-
  qc_summary_saved$MRPRESSO_global_p < 0.05


heterogeneity_flag_match <-
  identical(
    as.logical(
      qc_summary_saved$heterogeneity_evidence
    ),
    expected_heterogeneity
  )

egger_flag_match <-
  identical(
    as.logical(
      qc_summary_saved$directional_pleiotropy_evidence
    ),
    expected_egger_pleiotropy
  )

presso_flag_match <-
  identical(
    as.logical(
      qc_summary_saved$presso_global_evidence
    ),
    expected_presso
  )


cat(
  "Heterogeneity flags consistent: ",
  heterogeneity_flag_match,
  "\n"
)

cat(
  "Egger pleiotropy flags consistent: ",
  egger_flag_match,
  "\n"
)

cat(
  "MR-PRESSO flags consistent: ",
  presso_flag_match,
  "\n"
)
#确认最终evidence数
cat(
  "\nHeterogeneity evidence:",
  sum(expected_heterogeneity),
  "/ 10\n"
)

cat(
  "Directional pleiotropy evidence:",
  sum(expected_egger_pleiotropy),
  "/ 10\n"
)

cat(
  "MR-PRESSO global evidence:",
  sum(expected_presso),
  "/ 10\n"
)
# ============================================================
# 15H. Save QC provenance audit
# ============================================================

write.csv(
  source_qc,
  "00_admin/audit_outputs/source_QC_results.csv",
  row.names = FALSE
)

write.csv(
  qc_provenance_match,
  "00_admin/audit_outputs/QC_provenance_matching.csv",
  row.names = FALSE
)


cat("\nQC provenance outputs saved.\n")
# ============================================================
# 16. REPRODUCIBILITY AUDIT
# ============================================================

dir.create(
  "00_admin/audit_outputs",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 16A. Save complete R session information
# ============================================================

session_info_txt <- capture.output(
  sessionInfo()
)

writeLines(
  session_info_txt,
  "00_admin/audit_outputs/sessionInfo_2026-09-07.txt"
)

cat("\nR session information saved.\n")


# ============================================================
# 16B. R version
# ============================================================

r_version_audit <- data.frame(
  
  item = c(
    "R_version",
    "R_platform",
    "OS"
  ),
  
  value = c(
    R.version.string,
    R.version$platform,
    paste(
      Sys.info()[["sysname"]],
      Sys.info()[["release"]]
    )
  ),
  
  stringsAsFactors = FALSE
)

print(r_version_audit)

write.csv(
  r_version_audit,
  "00_admin/audit_outputs/R_environment.csv",
  row.names = FALSE
)


# ============================================================
# 16C. Identify packages referenced in project R scripts
# ============================================================

project_scripts <- list.files(
  "04_R_scripts",
  pattern = "\\.R$",
  full.names = TRUE
)

script_text <- unlist(
  lapply(
    project_scripts,
    readLines,
    warn = FALSE
  )
)


# library(package)
library_hits <- grep(
  "\\b(library|require)\\s*\\(",
  script_text,
  value = TRUE
)

library_packages <- gsub(
  ".*\\b(?:library|require)\\s*\\(\\s*[\"']?([^\"') ,]+).*",
  "\\1",
  library_hits,
  perl = TRUE
)


# package::function
namespace_hits <- regmatches(
  script_text,
  gregexpr(
    "[A-Za-z][A-Za-z0-9.]*::[A-Za-z][A-Za-z0-9._]*",
    script_text,
    perl = TRUE
  )
)

namespace_hits <- unlist(namespace_hits)

namespace_packages <- sub(
  "::.*$",
  "",
  namespace_hits
)


packages_used <- sort(
  unique(
    c(
      library_packages,
      namespace_packages
    )
  )
)

packages_used <- packages_used[
  nzchar(packages_used)
]

cat("\nPackages detected in R scripts:\n")
print(packages_used)


# ============================================================
# 16D. Record installed package versions
# ============================================================

package_versions <- data.frame(
  
  package = packages_used,
  
  installed = packages_used %in%
    rownames(installed.packages()),
  
  version = sapply(
    packages_used,
    function(pkg) {
      
      if (
        pkg %in% rownames(installed.packages())
      ) {
        
        as.character(
          packageVersion(pkg)
        )
        
      } else {
        
        NA_character_
      }
    }
  ),
  
  stringsAsFactors = FALSE
)


print(
  package_versions,
  row.names = FALSE
)

write.csv(
  package_versions,
  "00_admin/audit_outputs/package_versions.csv",
  row.names = FALSE
)


# ============================================================
# 16E. Create R-script manifest
# ============================================================

script_info <- file.info(
  project_scripts
)

script_manifest <- data.frame(
  
  script =
    basename(project_scripts),
  
  path =
    project_scripts,
  
  size_bytes =
    script_info$size,
  
  modified_time =
    script_info$mtime,
  
  md5 =
    unname(
      tools::md5sum(
        project_scripts
      )
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  script_manifest,
  "00_admin/audit_outputs/R_script_manifest.csv",
  row.names = FALSE
)


cat("\nR script manifest:\n")

print(
  script_manifest[
    ,
    c(
      "script",
      "md5"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 16F. Verify critical analysis scripts exist
# ============================================================

critical_scripts <- c(
  "00_functions.R",
  "04_primary_MR.R",
  "05_threshold_sensitivity.R",
  "06_finngen_sensitivity.R",
  "07_east_asian_sensitivity.R",
  "08_reverse_MR.R",
  "09_summary_results.R",
  "11_core_audit_checks.R"
)


critical_exists <- file.exists(
  file.path(
    "04_R_scripts",
    critical_scripts
  )
)


critical_script_check <- data.frame(
  
  script =
    critical_scripts,
  
  exists =
    critical_exists,
  
  stringsAsFactors = FALSE
)


cat("\nCritical script check:\n")

print(
  critical_script_check,
  row.names = FALSE
)

write.csv(
  critical_script_check,
  "00_admin/audit_outputs/critical_script_check.csv",
  row.names = FALSE
)


# ============================================================
# 16G. Verify critical frozen summary results
# ============================================================

critical_results <- c(
  
  "05_results/06_summary/all_MR_results.csv",
  
  "05_results/06_summary/IVW_summary.csv",
  
  "05_results/06_summary/IVW_summary_formatted.csv",
  
  "05_results/06_summary/QC_summary.csv"
)


critical_results_exist <-
  file.exists(
    critical_results
  )


critical_result_check <- data.frame(
  
  file =
    critical_results,
  
  exists =
    critical_results_exist,
  
  stringsAsFactors = FALSE
)


cat("\nCritical result-file check:\n")

print(
  critical_result_check,
  row.names = FALSE
)

write.csv(
  critical_result_check,
  "00_admin/audit_outputs/critical_result_check.csv",
  row.names = FALSE
)


# ============================================================
# 16H. Verify frozen snapshot
# ============================================================

snapshot_path <-
  "09_snapshot/2026-09-05_preregistered_core_v1.0"

snapshot_exists <-
  dir.exists(
    snapshot_path
  )

snapshot_md5_exists <-
  file.exists(
    file.path(
      snapshot_path,
      "snapshot_MD5_manifest.csv"
    )
  )

snapshot_readme_exists <-
  file.exists(
    file.path(
      snapshot_path,
      "README_snapshot.txt"
    )
  )


cat(
  "\nFrozen snapshot exists: ",
  snapshot_exists,
  "\n"
)

cat(
  "Snapshot MD5 manifest exists: ",
  snapshot_md5_exists,
  "\n"
)

cat(
  "Snapshot README exists: ",
  snapshot_readme_exists,
  "\n"
)


# ============================================================
# 16I. Final reproducibility checks
# ============================================================

all_critical_scripts_exist <-
  all(
    critical_script_check$exists
  )

all_critical_results_exist <-
  all(
    critical_result_check$exists
  )

all_detected_packages_installed <-
  all(
    package_versions$installed
  )


cat("\n========================================\n")
cat("REPRODUCIBILITY AUDIT CHECKS\n")
cat("========================================\n")

cat(
  "All critical R scripts exist: ",
  all_critical_scripts_exist,
  "\n"
)

cat(
  "All critical result files exist: ",
  all_critical_results_exist,
  "\n"
)

cat(
  "All detected R packages installed: ",
  all_detected_packages_installed,
  "\n"
)

cat(
  "Frozen core snapshot exists: ",
  snapshot_exists,
  "\n"
)

cat(
  "Snapshot MD5 manifest exists: ",
  snapshot_md5_exists,
  "\n"
)

cat(
  "Snapshot README exists: ",
  snapshot_readme_exists,
  "\n"
)