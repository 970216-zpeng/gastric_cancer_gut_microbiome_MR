# ============================================================
# 10_exploratory_exposure_replication.R
#
# Post-registration exploratory exposure-source replication
# Phase 1: Exposure-only feasibility audit
#
# OSF Secondary Data Preregistration:
# Kv7tj
#
# IMPORTANT:
# This phase does NOT access GCST90018849 outcome associations.
# ============================================================

#install.packages("R.utils")
# ============================================================
# 0. Project root
# ============================================================

if (basename(getwd()) == "04_R_scripts") {
  setwd("..")
}

cat("Working directory:\n")
print(getwd())


stopifnot(
  dir.exists("01_raw_data/01_exposure/Swedish2026"),
  dir.exists("01_raw_data/01_exposure/DMP2022")
)


# ============================================================
# 1. Packages
# ============================================================

library(data.table)
library(dplyr)


# ============================================================
# 2. Create exploratory intermediate-data folders
# ============================================================

dir.create(
  "02_intermediate_data/05_exposure_source_replication",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "02_intermediate_data/05_exposure_source_replication/Swedish2026",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "02_intermediate_data/05_exposure_source_replication/DMP2022",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Helper: locate one harmonised GWAS file by accession
# ============================================================

find_harmonised_gwas <- function(folder, accession) {
  
  x <- list.files(
    folder,
    full.names = TRUE
  )
  
  x <- x[
    grepl(
      accession,
      basename(x),
      fixed = TRUE
    )
  ]
  
  # keep harmonised summary-statistic file only
  # exclude YAML / metadata files
  x <- x[
    grepl(
      "\\.h\\.tsv(\\.gz)?$",
      basename(x),
      ignore.case = TRUE
    )
  ]
  
  if (length(x) != 1) {
    
    cat("\nFiles found for", accession, ":\n")
    print(x)
    
    stop(
      paste(
        "Expected exactly one harmonised GWAS file for",
        accession
      )
    )
  }
  
  normalizePath(
    x,
    winslash = "/",
    mustWork = TRUE
  )
}


# ============================================================
# 4. Locate the four registered exposure datasets
# ============================================================

swedish_v_file <- find_harmonised_gwas(
  "01_raw_data/01_exposure/Swedish2026",
  "GCST90671339"
)

swedish_va_file <- find_harmonised_gwas(
  "01_raw_data/01_exposure/Swedish2026",
  "GCST90671681"
)

dmp_v_file <- find_harmonised_gwas(
  "01_raw_data/01_exposure/DMP2022",
  "GCST90027724"
)

dmp_va_file <- find_harmonised_gwas(
  "01_raw_data/01_exposure/DMP2022",
  "GCST90027679"
)


cat("\n========================================\n")
cat("REGISTERED EXPOSURE FILES\n")
cat("========================================\n")

cat("\nSwedish Veillonella:\n", swedish_v_file, "\n")
cat("\nSwedish Veillonellaceae:\n", swedish_va_file, "\n")
cat("\nDMP Veillonella:\n", dmp_v_file, "\n")
cat("\nDMP Veillonellaceae:\n", dmp_va_file, "\n")


# ============================================================
# 5. File-size audit
# ============================================================

exposure_file_manifest <- data.frame(
  
  dataset = c(
    "Swedish2026",
    "Swedish2026",
    "DMP2022",
    "DMP2022"
  ),
  
  taxon = c(
    "Veillonella",
    "Veillonellaceae",
    "Veillonella",
    "Veillonellaceae"
  ),
  
  accession = c(
    "GCST90671339",
    "GCST90671681",
    "GCST90027724",
    "GCST90027679"
  ),
  
  file = c(
    swedish_v_file,
    swedish_va_file,
    dmp_v_file,
    dmp_va_file
  ),
  
  stringsAsFactors = FALSE
)


exposure_file_manifest$size_MB <-
  round(
    file.info(
      exposure_file_manifest$file
    )$size / 1024^2,
    1
  )


print(
  exposure_file_manifest[
    ,
    c(
      "dataset",
      "taxon",
      "accession",
      "size_MB"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 6. Helper: inspect harmonised GWAS header only
# ============================================================

inspect_gwas_header <- function(path, label) {
  
  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n")
  
  x <- data.table::fread(
    path,
    nrows = 5,
    showProgress = FALSE
  )
  
  cat("\nColumn names:\n")
  print(names(x))
  
  cat("\nFirst rows:\n")
  print(x)
  
  invisible(x)
}


# ============================================================
# 7. Inspect all four files
# ============================================================

head_swedish_v <- inspect_gwas_header(
  swedish_v_file,
  "Swedish 2026 - Veillonella"
)

head_swedish_va <- inspect_gwas_header(
  swedish_va_file,
  "Swedish 2026 - Veillonellaceae"
)

head_dmp_v <- inspect_gwas_header(
  dmp_v_file,
  "DMP 2022 - Veillonella"
)

head_dmp_va <- inspect_gwas_header(
  dmp_va_file,
  "DMP 2022 - Veillonellaceae"
)


# ============================================================
# 8. Save file manifest
# ============================================================

write.csv(
  exposure_file_manifest,
  "02_intermediate_data/05_exposure_source_replication/exposure_file_manifest.csv",
  row.names = FALSE
)

cat(
  "\nSaved:\n",
  "02_intermediate_data/05_exposure_source_replication/exposure_file_manifest.csv\n"
)
# ============================================================
# 9. Compact column-name audit
# ============================================================

cat("\n===== Swedish 2026: Veillonella =====\n")
print(names(head_swedish_v))

cat("\n===== Swedish 2026: Veillonellaceae =====\n")
print(names(head_swedish_va))

cat("\n===== DMP 2022: Veillonella =====\n")
print(names(head_dmp_v))

cat("\n===== DMP 2022: Veillonellaceae =====\n")
print(names(head_dmp_va))
#检查
identical(
  names(head_swedish_v),
  names(head_swedish_va)
)

identical(
  names(head_swedish_v),
  names(head_dmp_v)
)

identical(
  names(head_dmp_v),
  names(head_dmp_va)
)
# ============================================================
# 10. Helper: identify required columns automatically
# ============================================================

pick_col <- function(nms, candidates, label) {
  
  hit <- candidates[candidates %in% nms]
  
  if (length(hit) == 0) {
    stop(
      paste(
        "Could not identify column for:",
        label
      )
    )
  }
  
  hit[1]
}


get_gwas_mapping <- function(header_object) {
  
  nms <- names(header_object)
  
  data.frame(
    
    field = c(
      "SNP",
      "beta",
      "se",
      "pval",
      "effect_allele",
      "other_allele",
      "eaf",
      "chr",
      "pos"
    ),
    
    column = c(
      
      pick_col(
        nms,
        c(
          "hm_rsid",
          "rsid",
          "rs_id",
          "SNP",
          "variant_id"
        ),
        "SNP / rsID"
      ),
      
      pick_col(
        nms,
        c(
          "beta",
          "hm_beta",
          "effect_size",
          "BETA"
        ),
        "beta"
      ),
      
      pick_col(
        nms,
        c(
          "standard_error",
          "se",
          "hm_standard_error",
          "SE"
        ),
        "standard error"
      ),
      
      pick_col(
        nms,
        c(
          "p_value",
          "pval",
          "p_value",
          "P"
        ),
        "P value"
      ),
      
      pick_col(
        nms,
        c(
          "effect_allele",
          "hm_effect_allele"
        ),
        "effect allele"
      ),
      
      pick_col(
        nms,
        c(
          "other_allele",
          "hm_other_allele"
        ),
        "other allele"
      ),
      
      pick_col(
        nms,
        c(
          "effect_allele_frequency",
          "hm_effect_allele_frequency",
          "eaf"
        ),
        "effect allele frequency"
      ),
      
      pick_col(
        nms,
        c(
          "chromosome",
          "hm_chrom",
          "chr"
        ),
        "chromosome"
      ),
      
      pick_col(
        nms,
        c(
          "base_pair_location",
          "hm_pos",
          "position",
          "pos"
        ),
        "position"
      )
    ),
    
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 11. Build column mappings
# ============================================================

map_swedish <- get_gwas_mapping(
  head_swedish_v
)

map_dmp <- get_gwas_mapping(
  head_dmp_v
)

cat("\n========================================\n")
cat("SWEDISH COLUMN MAPPING\n")
cat("========================================\n")

print(
  map_swedish,
  row.names = FALSE
)


cat("\n========================================\n")
cat("DMP COLUMN MAPPING\n")
cat("========================================\n")

print(
  map_dmp,
  row.names = FALSE
)
# ============================================================
# 12. Helper: extract candidate IVs from one exposure GWAS
# ============================================================

extract_candidate_ivs <- function(
    path,
    mapping,
    p_threshold,
    dataset,
    taxon,
    accession
) {
  
  selected_cols <- unique(
    mapping$column
  )
  
  cat(
    "\nReading:",
    dataset,
    taxon,
    accession,
    "\n"
  )
  
  x <- data.table::fread(
    path,
    select = selected_cols,
    showProgress = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Rename into common analytical schema
  # ----------------------------------------------------------
  
  get_col <- function(field) {
    mapping$column[
      mapping$field == field
    ]
  }
  
  
  out <- data.frame(
    
    SNP =
      as.character(
        x[[get_col("SNP")]]
      ),
    
    beta =
      as.numeric(
        x[[get_col("beta")]]
      ),
    
    se =
      as.numeric(
        x[[get_col("se")]]
      ),
    
    pval =
      as.numeric(
        x[[get_col("pval")]]
      ),
    
    effect_allele =
      as.character(
        x[[get_col("effect_allele")]]
      ),
    
    other_allele =
      as.character(
        x[[get_col("other_allele")]]
      ),
    
    eaf =
      as.numeric(
        x[[get_col("eaf")]]
      ),
    
    chr =
      as.character(
        x[[get_col("chr")]]
      ),
    
    pos =
      as.numeric(
        x[[get_col("pos")]]
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Initial QC
  # ----------------------------------------------------------
  
  n_total <- nrow(out)
  
  out <- out[
    !is.na(out$pval) &
      out$pval <= p_threshold,
    ,
    drop = FALSE
  ]
  
  n_p <- nrow(out)
  
  
  # Require essential MR fields
  out <- out[
    !is.na(out$SNP) &
      out$SNP != "" &
      grepl("^rs[0-9]+$", out$SNP) &
      !is.na(out$beta) &
      !is.na(out$se) &
      out$se > 0 &
      !is.na(out$effect_allele) &
      !is.na(out$other_allele),
    ,
    drop = FALSE
  ]
  
  n_complete <- nrow(out)
  
  
  # ----------------------------------------------------------
  # Resolve duplicated rsIDs
  # Keep strongest association
  # ----------------------------------------------------------
  
  out <- out[
    order(
      out$pval,
      decreasing = FALSE
    ),
    ,
    drop = FALSE
  ]
  
  out <- out[
    !duplicated(out$SNP),
    ,
    drop = FALSE
  ]
  
  n_unique <- nrow(out)
  
  
  # ----------------------------------------------------------
  # F statistic
  # ----------------------------------------------------------
  
  out$F_stat <-
    (out$beta^2) /
    (out$se^2)
  
  
  out <- out[
    !is.na(out$F_stat) &
      is.finite(out$F_stat) &
      out$F_stat > 10,
    ,
    drop = FALSE
  ]
  
  n_F <- nrow(out)
  
  
  out$dataset <- dataset
  out$taxon <- taxon
  out$accession <- accession
  out$p_threshold <- p_threshold
  
  
  audit <- data.frame(
    
    dataset = dataset,
    taxon = taxon,
    accession = accession,
    p_threshold = p_threshold,
    
    total_rows = n_total,
    p_threshold_rows = n_p,
    complete_valid_rows = n_complete,
    unique_rsids = n_unique,
    F_gt_10_rows = n_F,
    
    min_F_after_filter =
      if (n_F > 0)
        min(out$F_stat)
    else
      NA_real_,
    
    stringsAsFactors = FALSE
  )
  
  
  list(
    data = out,
    audit = audit
  )
}
# ============================================================
# 13. Primary exposure candidate IVs
#     P <= 5e-6
# ============================================================

swedish_v_p5e6 <- extract_candidate_ivs(
  path = swedish_v_file,
  mapping = map_swedish,
  p_threshold = 5e-6,
  dataset = "Swedish2026",
  taxon = "Veillonella",
  accession = "GCST90671339"
)

swedish_va_p5e6 <- extract_candidate_ivs(
  path = swedish_va_file,
  mapping = map_swedish,
  p_threshold = 5e-6,
  dataset = "Swedish2026",
  taxon = "Veillonellaceae",
  accession = "GCST90671681"
)

dmp_v_p5e6 <- extract_candidate_ivs(
  path = dmp_v_file,
  mapping = map_dmp,
  p_threshold = 5e-6,
  dataset = "DMP2022",
  taxon = "Veillonella",
  accession = "GCST90027724"
)

dmp_va_p5e6 <- extract_candidate_ivs(
  path = dmp_va_file,
  mapping = map_dmp,
  p_threshold = 5e-6,
  dataset = "DMP2022",
  taxon = "Veillonellaceae",
  accession = "GCST90027679"
)
# 14汇总
primary_candidate_audit <- dplyr::bind_rows(
  
  swedish_v_p5e6$audit,
  swedish_va_p5e6$audit,
  dmp_v_p5e6$audit,
  dmp_va_p5e6$audit
)


cat("\n========================================\n")
cat("PRIMARY P <= 5e-6 EXPOSURE FEASIBILITY\n")
cat("========================================\n")

print(
  primary_candidate_audit,
  row.names = FALSE
)
# 15保存
saveRDS(
  swedish_v_p5e6$data,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonella_candidate_p5e6_preclump.rds"
)

saveRDS(
  swedish_va_p5e6$data,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonellaceae_candidate_p5e6_preclump.rds"
)

saveRDS(
  dmp_v_p5e6$data,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonella_candidate_p5e6_preclump.rds"
)

saveRDS(
  dmp_va_p5e6$data,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonellaceae_candidate_p5e6_preclump.rds"
)


write.csv(
  primary_candidate_audit,
  "02_intermediate_data/05_exposure_source_replication/primary_p5e6_preclump_audit.csv",
  row.names = FALSE
)
# ============================================================
# 16. Prepare exposure data for LD clumping
# ============================================================

prepare_for_clumping <- function(x, exposure_name) {
  
  out <- data.frame(
    SNP = x$SNP,
    pval.exposure = x$pval,
    id.exposure = exposure_name,
    stringsAsFactors = FALSE
  )
  
  out
}


swedish_v_clump_input <- prepare_for_clumping(
  swedish_v_p5e6$data,
  "Swedish2026_Veillonella"
)

swedish_va_clump_input <- prepare_for_clumping(
  swedish_va_p5e6$data,
  "Swedish2026_Veillonellaceae"
)

dmp_v_clump_input <- prepare_for_clumping(
  dmp_v_p5e6$data,
  "DMP2022_Veillonella"
)

dmp_va_clump_input <- prepare_for_clumping(
  dmp_va_p5e6$data,
  "DMP2022_Veillonellaceae"
)
# ============================================================
# 17. LD clumping
#     Registered rule:
#     r2 < 0.001
#     10,000 kb
#     European LD reference
# ============================================================

swedish_v_clumped <- TwoSampleMR::clump_data(
  swedish_v_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-6,
  pop = "EUR"
)

swedish_va_clumped <- TwoSampleMR::clump_data(
  swedish_va_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-6,
  pop = "EUR"
)

dmp_v_clumped <- TwoSampleMR::clump_data(
  dmp_v_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-6,
  pop = "EUR"
)

dmp_va_clumped <- TwoSampleMR::clump_data(
  dmp_va_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-6,
  pop = "EUR"
)
# ============================================================
# 18. Recover complete exposure information after clumping
# ============================================================

recover_clumped_data <- function(candidate_data, clumped_data) {
  
  out <- candidate_data[
    candidate_data$SNP %in% clumped_data$SNP,
    ,
    drop = FALSE
  ]
  
  out <- out[
    match(
      clumped_data$SNP,
      out$SNP
    ),
    ,
    drop = FALSE
  ]
  
  rownames(out) <- NULL
  
  out
}


swedish_v_primary_IV <- recover_clumped_data(
  swedish_v_p5e6$data,
  swedish_v_clumped
)

swedish_va_primary_IV <- recover_clumped_data(
  swedish_va_p5e6$data,
  swedish_va_clumped
)

dmp_v_primary_IV <- recover_clumped_data(
  dmp_v_p5e6$data,
  dmp_v_clumped
)

dmp_va_primary_IV <- recover_clumped_data(
  dmp_va_p5e6$data,
  dmp_va_clumped
)
# ============================================================
# 19. Final primary IV counts after LD clumping
# ============================================================

primary_IV_audit <- data.frame(
  
  dataset = c(
    "Swedish2026",
    "Swedish2026",
    "DMP2022",
    "DMP2022"
  ),
  
  taxon = c(
    "Veillonella",
    "Veillonellaceae",
    "Veillonella",
    "Veillonellaceae"
  ),
  
  preclump_F_gt_10 = c(
    nrow(swedish_v_p5e6$data),
    nrow(swedish_va_p5e6$data),
    nrow(dmp_v_p5e6$data),
    nrow(dmp_va_p5e6$data)
  ),
  
  final_clumped_IVs = c(
    nrow(swedish_v_primary_IV),
    nrow(swedish_va_primary_IV),
    nrow(dmp_v_primary_IV),
    nrow(dmp_va_primary_IV)
  ),
  
  min_F = c(
    min(swedish_v_primary_IV$F_stat),
    min(swedish_va_primary_IV$F_stat),
    min(dmp_v_primary_IV$F_stat),
    min(dmp_va_primary_IV$F_stat)
  ),
  
  median_F = c(
    median(swedish_v_primary_IV$F_stat),
    median(swedish_va_primary_IV$F_stat),
    median(dmp_v_primary_IV$F_stat),
    median(dmp_va_primary_IV$F_stat)
  ),
  
  max_F = c(
    max(swedish_v_primary_IV$F_stat),
    max(swedish_va_primary_IV$F_stat),
    max(dmp_v_primary_IV$F_stat),
    max(dmp_va_primary_IV$F_stat)
  ),
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("PRIMARY EXTERNAL EXPOSURE IV SETS\n")
cat("========================================\n")

print(
  primary_IV_audit,
  row.names = FALSE
)
# ============================================================
# 20. Save primary clumped IV sets
# ============================================================

saveRDS(
  swedish_v_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonella_IV_primary_p5e6.rds"
)

write.csv(
  swedish_v_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonella_IV_primary_p5e6.csv",
  row.names = FALSE
)


saveRDS(
  swedish_va_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonellaceae_IV_primary_p5e6.rds"
)

write.csv(
  swedish_va_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonellaceae_IV_primary_p5e6.csv",
  row.names = FALSE
)


saveRDS(
  dmp_v_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonella_IV_primary_p5e6.rds"
)

write.csv(
  dmp_v_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonella_IV_primary_p5e6.csv",
  row.names = FALSE
)


saveRDS(
  dmp_va_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonellaceae_IV_primary_p5e6.rds"
)

write.csv(
  dmp_va_primary_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonellaceae_IV_primary_p5e6.csv",
  row.names = FALSE
)


write.csv(
  primary_IV_audit,
  "02_intermediate_data/05_exposure_source_replication/primary_p5e6_final_IV_audit.csv",
  row.names = FALSE
)
# ============================================================
# 21. Threshold-sensitivity candidate IVs
#     P <= 1e-5
# ============================================================

swedish_v_p1e5 <- extract_candidate_ivs(
  path = swedish_v_file,
  mapping = map_swedish,
  p_threshold = 1e-5,
  dataset = "Swedish2026",
  taxon = "Veillonella",
  accession = "GCST90671339"
)

swedish_va_p1e5 <- extract_candidate_ivs(
  path = swedish_va_file,
  mapping = map_swedish,
  p_threshold = 1e-5,
  dataset = "Swedish2026",
  taxon = "Veillonellaceae",
  accession = "GCST90671681"
)

dmp_v_p1e5 <- extract_candidate_ivs(
  path = dmp_v_file,
  mapping = map_dmp,
  p_threshold = 1e-5,
  dataset = "DMP2022",
  taxon = "Veillonella",
  accession = "GCST90027724"
)

dmp_va_p1e5 <- extract_candidate_ivs(
  path = dmp_va_file,
  mapping = map_dmp,
  p_threshold = 1e-5,
  dataset = "DMP2022",
  taxon = "Veillonellaceae",
  accession = "GCST90027679"
)
# ============================================================
# 22. Threshold pre-clump feasibility audit
# ============================================================

threshold_candidate_audit <- dplyr::bind_rows(
  swedish_v_p1e5$audit,
  swedish_va_p1e5$audit,
  dmp_v_p1e5$audit,
  dmp_va_p1e5$audit
)

cat("\n========================================\n")
cat("THRESHOLD P <= 1e-5 PRECLUMP FEASIBILITY\n")
cat("========================================\n")

print(
  threshold_candidate_audit,
  row.names = FALSE
)
# ============================================================
# 23. Prepare threshold-sensitivity clumping input
# ============================================================

swedish_v_p1e5_clump_input <- prepare_for_clumping(
  swedish_v_p1e5$data,
  "Swedish2026_Veillonella_p1e5"
)

swedish_va_p1e5_clump_input <- prepare_for_clumping(
  swedish_va_p1e5$data,
  "Swedish2026_Veillonellaceae_p1e5"
)

dmp_v_p1e5_clump_input <- prepare_for_clumping(
  dmp_v_p1e5$data,
  "DMP2022_Veillonella_p1e5"
)

dmp_va_p1e5_clump_input <- prepare_for_clumping(
  dmp_va_p1e5$data,
  "DMP2022_Veillonellaceae_p1e5"
)
# ============================================================
# 24. Threshold-sensitivity LD clumping
#
# P <= 1e-5
# r2 < 0.001
# 10,000 kb
# EUR reference
# ============================================================

swedish_v_p1e5_clumped <- TwoSampleMR::clump_data(
  swedish_v_p1e5_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-5,
  pop = "EUR"
)

swedish_va_p1e5_clumped <- TwoSampleMR::clump_data(
  swedish_va_p1e5_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-5,
  pop = "EUR"
)

dmp_v_p1e5_clumped <- TwoSampleMR::clump_data(
  dmp_v_p1e5_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-5,
  pop = "EUR"
)

dmp_va_p1e5_clumped <- TwoSampleMR::clump_data(
  dmp_va_p1e5_clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-5,
  pop = "EUR"
)
# ============================================================
# 25. Recover complete threshold-sensitivity IV information
# ============================================================

swedish_v_sensitivity_IV <- recover_clumped_data(
  swedish_v_p1e5$data,
  swedish_v_p1e5_clumped
)

swedish_va_sensitivity_IV <- recover_clumped_data(
  swedish_va_p1e5$data,
  swedish_va_p1e5_clumped
)

dmp_v_sensitivity_IV <- recover_clumped_data(
  dmp_v_p1e5$data,
  dmp_v_p1e5_clumped
)

dmp_va_sensitivity_IV <- recover_clumped_data(
  dmp_va_p1e5$data,
  dmp_va_p1e5_clumped
)
# ============================================================
# 26. Final threshold-sensitivity IV audit
# ============================================================

threshold_IV_audit <- data.frame(
  
  dataset = c(
    "Swedish2026",
    "Swedish2026",
    "DMP2022",
    "DMP2022"
  ),
  
  taxon = c(
    "Veillonella",
    "Veillonellaceae",
    "Veillonella",
    "Veillonellaceae"
  ),
  
  preclump_F_gt_10 = c(
    nrow(swedish_v_p1e5$data),
    nrow(swedish_va_p1e5$data),
    nrow(dmp_v_p1e5$data),
    nrow(dmp_va_p1e5$data)
  ),
  
  final_clumped_IVs = c(
    nrow(swedish_v_sensitivity_IV),
    nrow(swedish_va_sensitivity_IV),
    nrow(dmp_v_sensitivity_IV),
    nrow(dmp_va_sensitivity_IV)
  ),
  
  min_F = c(
    min(swedish_v_sensitivity_IV$F_stat),
    min(swedish_va_sensitivity_IV$F_stat),
    min(dmp_v_sensitivity_IV$F_stat),
    min(dmp_va_sensitivity_IV$F_stat)
  ),
  
  median_F = c(
    median(swedish_v_sensitivity_IV$F_stat),
    median(swedish_va_sensitivity_IV$F_stat),
    median(dmp_v_sensitivity_IV$F_stat),
    median(dmp_va_sensitivity_IV$F_stat)
  ),
  
  max_F = c(
    max(swedish_v_sensitivity_IV$F_stat),
    max(swedish_va_sensitivity_IV$F_stat),
    max(dmp_v_sensitivity_IV$F_stat),
    max(dmp_va_sensitivity_IV$F_stat)
  ),
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("THRESHOLD P <= 1e-5 EXTERNAL EXPOSURE IV SETS\n")
cat("========================================\n")

print(
  threshold_IV_audit,
  row.names = FALSE
)
# ============================================================
# 27. Save threshold-sensitivity IV sets
# ============================================================

saveRDS(
  swedish_v_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonella_IV_sensitivity_p1e5.rds"
)

write.csv(
  swedish_v_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonella_IV_sensitivity_p1e5.csv",
  row.names = FALSE
)


saveRDS(
  swedish_va_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonellaceae_IV_sensitivity_p1e5.rds"
)

write.csv(
  swedish_va_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/Swedish2026/Veillonellaceae_IV_sensitivity_p1e5.csv",
  row.names = FALSE
)


saveRDS(
  dmp_v_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonella_IV_sensitivity_p1e5.rds"
)

write.csv(
  dmp_v_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonella_IV_sensitivity_p1e5.csv",
  row.names = FALSE
)


saveRDS(
  dmp_va_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonellaceae_IV_sensitivity_p1e5.rds"
)

write.csv(
  dmp_va_sensitivity_IV,
  "02_intermediate_data/05_exposure_source_replication/DMP2022/Veillonellaceae_IV_sensitivity_p1e5.csv",
  row.names = FALSE
)


write.csv(
  threshold_IV_audit,
  "02_intermediate_data/05_exposure_source_replication/threshold_p1e5_final_IV_audit.csv",
  row.names = FALSE
)
# ============================================================
# 28. Build union of all registered external exposure IVs
# ============================================================

all_external_target_snps <- unique(
  c(
    swedish_v_primary_IV$SNP,
    swedish_va_primary_IV$SNP,
    dmp_v_primary_IV$SNP,
    dmp_va_primary_IV$SNP,
    
    swedish_v_sensitivity_IV$SNP,
    swedish_va_sensitivity_IV$SNP,
    dmp_v_sensitivity_IV$SNP,
    dmp_va_sensitivity_IV$SNP
  )
)

all_external_target_snps <- all_external_target_snps[
  !is.na(all_external_target_snps) &
    all_external_target_snps != ""
]

cat("\n========================================\n")
cat("EXTERNAL EXPOSURE TARGET SNP UNION\n")
cat("========================================\n")

cat(
  "Unique SNPs to extract from GCST90018849:",
  length(all_external_target_snps),
  "\n"
)

print(all_external_target_snps)
# ============================================================
# 29. Load reusable project functions
# ============================================================

source(
  "04_R_scripts/00_functions.R"
)

stopifnot(
  exists("extract_outcome_snps")
)
# ============================================================
# 30. Primary gastric cancer outcome
# ============================================================

gc_outcome_file <-
  "01_raw_data/02_outcome/ebi-a-GCST90018849.vcf.gz"

stopifnot(
  file.exists(gc_outcome_file)
)

cat(
  "\nOutcome file:\n",
  gc_outcome_file,
  "\n"
)
# ============================================================
# 31. Extract GCST90018849 associations
#     for all registered external exposure IVs
# ============================================================

external_gc_extract <- extract_outcome_snps(
  vcf_file = gc_outcome_file,
  target_snps = all_external_target_snps,
  chunk_size = 50000
)
# ============================================================
# 32. Outcome extraction availability audit
# ============================================================

cat("\n========================================\n")
cat("GCST90018849 OUTCOME EXTRACTION AUDIT\n")
cat("========================================\n")

cat(
  "Requested:",
  external_gc_extract$n_requested,
  "\n"
)

cat(
  "Found:",
  external_gc_extract$n_found,
  "\n"
)

cat(
  "Missing:",
  external_gc_extract$n_missing,
  "\n"
)

cat("\nMissing SNPs:\n")
print(
  external_gc_extract$missing_snps
)
#提取outcome数据本身
external_gc_raw <- external_gc_extract$data

cat("\nOutcome rows:\n")
print(nrow(external_gc_raw))

cat("\nColumn names:\n")
print(names(external_gc_raw))

cat("\nMissing beta:\n")
print(sum(is.na(external_gc_raw$beta)))

cat("\nMissing se:\n")
print(sum(is.na(external_gc_raw$se)))

cat("\nMissing pval:\n")
print(sum(is.na(external_gc_raw$pval)))

cat("\nDuplicate SNPs:\n")
print(sum(duplicated(external_gc_raw$SNP)))
# ============================================================
# 33. Save extracted external-replication outcome data
# ============================================================

dir.create(
  "02_intermediate_data/05_exposure_source_replication/GCST90018849",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  external_gc_raw,
  "02_intermediate_data/05_exposure_source_replication/GCST90018849/all_external_IVs_outcome_extracted.rds"
)

write.csv(
  external_gc_raw,
  "02_intermediate_data/05_exposure_source_replication/GCST90018849/all_external_IVs_outcome_extracted.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    SNP = external_gc_extract$missing_snps
  ),
  "02_intermediate_data/05_exposure_source_replication/GCST90018849/all_external_IVs_missing_in_outcome.csv",
  row.names = FALSE
)
# ============================================================
# 34. Analysis-specific outcome availability audit
# ============================================================

audit_outcome_availability <- function(iv_data, label) {
  
  available_snps <- intersect(
    iv_data$SNP,
    external_gc_raw$SNP
  )
  
  missing_snps <- setdiff(
    iv_data$SNP,
    external_gc_raw$SNP
  )
  
  data.frame(
    analysis = label,
    selected_IVs = nrow(iv_data),
    outcome_available = length(available_snps),
    outcome_missing = length(missing_snps),
    stringsAsFactors = FALSE
  )
}


external_outcome_flow <- dplyr::bind_rows(
  
  audit_outcome_availability(
    swedish_v_primary_IV,
    "Swedish2026_Veillonella_primary_p5e6"
  ),
  
  audit_outcome_availability(
    swedish_va_primary_IV,
    "Swedish2026_Veillonellaceae_primary_p5e6"
  ),
  
  audit_outcome_availability(
    dmp_v_primary_IV,
    "DMP2022_Veillonella_primary_p5e6"
  ),
  
  audit_outcome_availability(
    dmp_va_primary_IV,
    "DMP2022_Veillonellaceae_primary_p5e6"
  ),
  
  audit_outcome_availability(
    swedish_v_sensitivity_IV,
    "Swedish2026_Veillonella_sensitivity_p1e5"
  ),
  
  audit_outcome_availability(
    swedish_va_sensitivity_IV,
    "Swedish2026_Veillonellaceae_sensitivity_p1e5"
  ),
  
  audit_outcome_availability(
    dmp_v_sensitivity_IV,
    "DMP2022_Veillonella_sensitivity_p1e5"
  ),
  
  audit_outcome_availability(
    dmp_va_sensitivity_IV,
    "DMP2022_Veillonellaceae_sensitivity_p1e5"
  )
)


cat("\n========================================\n")
cat("ANALYSIS-SPECIFIC OUTCOME AVAILABILITY\n")
cat("========================================\n")

print(
  external_outcome_flow,
  row.names = FALSE
)


write.csv(
  external_outcome_flow,
  "02_intermediate_data/05_exposure_source_replication/GCST90018849/analysis_specific_outcome_availability.csv",
  row.names = FALSE
)
# ============================================================
# 35. Helper functions for harmonisation
# ============================================================

make_exposure_dat <- function(iv_data, exposure_label, exposure_id) {
  
  data.frame(
    
    SNP = iv_data$SNP,
    
    beta.exposure = iv_data$beta,
    se.exposure = iv_data$se,
    effect_allele.exposure = iv_data$effect_allele,
    other_allele.exposure = iv_data$other_allele,
    eaf.exposure = iv_data$eaf,
    pval.exposure = iv_data$pval,
    
    exposure = exposure_label,
    id.exposure = exposure_id,
    
    stringsAsFactors = FALSE
  )
}


make_outcome_dat <- function(target_snps) {
  
  x <- external_gc_raw[
    external_gc_raw$SNP %in% target_snps,
    ,
    drop = FALSE
  ]
  
  data.frame(
    
    SNP = x$SNP,
    
    beta.outcome = x$beta,
    se.outcome = x$se,
    effect_allele.outcome = x$effect_allele,
    other_allele.outcome = x$other_allele,
    eaf.outcome = x$eaf,
    pval.outcome = x$pval,
    
    outcome = "Gastric cancer",
    id.outcome = "GCST90018849",
    
    samplesize.outcome = x$samplesize,
    ncase.outcome = x$ncase,
    
    stringsAsFactors = FALSE
  )
}
# ============================================================
# 36. Harmonisation helper
# ============================================================

harmonise_external_analysis <- function(
    iv_data,
    exposure_label,
    exposure_id
) {
  
  exp_dat <- make_exposure_dat(
    iv_data,
    exposure_label,
    exposure_id
  )
  
  out_dat <- make_outcome_dat(
    iv_data$SNP
  )
  
  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat,
    action = 3
  )
  
  harmonised
}
# ============================================================
# 37. Primary P <= 5e-6 harmonisation
# ============================================================

harm_swedish_v_primary <- harmonise_external_analysis(
  swedish_v_primary_IV,
  "Swedish2026 Veillonella",
  "GCST90671339"
)

harm_swedish_va_primary <- harmonise_external_analysis(
  swedish_va_primary_IV,
  "Swedish2026 Veillonellaceae",
  "GCST90671681"
)

harm_dmp_v_primary <- harmonise_external_analysis(
  dmp_v_primary_IV,
  "DMP2022 Veillonella",
  "GCST90027724"
)

harm_dmp_va_primary <- harmonise_external_analysis(
  dmp_va_primary_IV,
  "DMP2022 Veillonellaceae",
  "GCST90027679"
)
# ============================================================
# 38. Sensitivity P <= 1e-5 harmonisation
# ============================================================

harm_swedish_v_p1e5 <- harmonise_external_analysis(
  swedish_v_sensitivity_IV,
  "Swedish2026 Veillonella p1e-5",
  "GCST90671339"
)

harm_swedish_va_p1e5 <- harmonise_external_analysis(
  swedish_va_sensitivity_IV,
  "Swedish2026 Veillonellaceae p1e-5",
  "GCST90671681"
)

harm_dmp_v_p1e5 <- harmonise_external_analysis(
  dmp_v_sensitivity_IV,
  "DMP2022 Veillonella p1e-5",
  "GCST90027724"
)

harm_dmp_va_p1e5 <- harmonise_external_analysis(
  dmp_va_sensitivity_IV,
  "DMP2022 Veillonellaceae p1e-5",
  "GCST90027679"
)
# ============================================================
# 39. Harmonisation audit
# ============================================================

audit_harmonisation <- function(dat, analysis_label) {
  
  data.frame(
    
    analysis = analysis_label,
    
    before_harmonisation =
      nrow(dat),
    
    retained_for_MR =
      sum(
        dat$mr_keep,
        na.rm = TRUE
      ),
    
    excluded =
      sum(
        !dat$mr_keep,
        na.rm = TRUE
      ),
    
    palindromic =
      if ("palindromic" %in% names(dat))
        sum(dat$palindromic, na.rm = TRUE)
    else
      NA_integer_,
    
    ambiguous =
      if ("ambiguous" %in% names(dat))
        sum(dat$ambiguous, na.rm = TRUE)
    else
      NA_integer_,
    
    stringsAsFactors = FALSE
  )
}


external_harmonisation_audit <- dplyr::bind_rows(
  
  audit_harmonisation(
    harm_swedish_v_primary,
    "Swedish2026_Veillonella_primary_p5e6"
  ),
  
  audit_harmonisation(
    harm_swedish_va_primary,
    "Swedish2026_Veillonellaceae_primary_p5e6"
  ),
  
  audit_harmonisation(
    harm_dmp_v_primary,
    "DMP2022_Veillonella_primary_p5e6"
  ),
  
  audit_harmonisation(
    harm_dmp_va_primary,
    "DMP2022_Veillonellaceae_primary_p5e6"
  ),
  
  audit_harmonisation(
    harm_swedish_v_p1e5,
    "Swedish2026_Veillonella_sensitivity_p1e5"
  ),
  
  audit_harmonisation(
    harm_swedish_va_p1e5,
    "Swedish2026_Veillonellaceae_sensitivity_p1e5"
  ),
  
  audit_harmonisation(
    harm_dmp_v_p1e5,
    "DMP2022_Veillonella_sensitivity_p1e5"
  ),
  
  audit_harmonisation(
    harm_dmp_va_p1e5,
    "DMP2022_Veillonellaceae_sensitivity_p1e5"
  )
)


cat("\n========================================\n")
cat("EXTERNAL EXPOSURE HARMONISATION AUDIT\n")
cat("========================================\n")

print(
  external_harmonisation_audit,
  row.names = FALSE
)
# ============================================================
# 40. MR-ready datasets
# ============================================================

mr_swedish_v_primary <- harm_swedish_v_primary[
  harm_swedish_v_primary$mr_keep,
]

mr_swedish_va_primary <- harm_swedish_va_primary[
  harm_swedish_va_primary$mr_keep,
]

mr_dmp_v_primary <- harm_dmp_v_primary[
  harm_dmp_v_primary$mr_keep,
]

mr_dmp_va_primary <- harm_dmp_va_primary[
  harm_dmp_va_primary$mr_keep,
]


mr_swedish_v_p1e5 <- harm_swedish_v_p1e5[
  harm_swedish_v_p1e5$mr_keep,
]

mr_swedish_va_p1e5 <- harm_swedish_va_p1e5[
  harm_swedish_va_p1e5$mr_keep,
]

mr_dmp_v_p1e5 <- harm_dmp_v_p1e5[
  harm_dmp_v_p1e5$mr_keep,
]

mr_dmp_va_p1e5 <- harm_dmp_va_p1e5[
  harm_dmp_va_p1e5$mr_keep,
]
#检查
cat("\n========================================\n")
cat("FINAL MR-READY SNP COUNTS\n")
cat("========================================\n")

cat("Swedish Veillonella primary:      ", nrow(mr_swedish_v_primary), "\n")
cat("Swedish Veillonellaceae primary:  ", nrow(mr_swedish_va_primary), "\n")
cat("DMP Veillonella primary:          ", nrow(mr_dmp_v_primary), "\n")
cat("DMP Veillonellaceae primary:      ", nrow(mr_dmp_va_primary), "\n")

cat("Swedish Veillonella p1e-5:        ", nrow(mr_swedish_v_p1e5), "\n")
cat("Swedish Veillonellaceae p1e-5:    ", nrow(mr_swedish_va_p1e5), "\n")
cat("DMP Veillonella p1e-5:            ", nrow(mr_dmp_v_p1e5), "\n")
cat("DMP Veillonellaceae p1e-5:        ", nrow(mr_dmp_va_p1e5), "\n")
# ============================================================
# 41. Save harmonisation audit and MR-ready datasets
# ============================================================

write.csv(
  external_harmonisation_audit,
  "02_intermediate_data/05_exposure_source_replication/GCST90018849/external_harmonisation_audit.csv",
  row.names = FALSE
)

save_external_harm <- function(dat, name) {
  
  saveRDS(
    dat,
    paste0(
      "02_intermediate_data/05_exposure_source_replication/GCST90018849/",
      name,
      ".rds"
    )
  )
  
  write.csv(
    dat,
    paste0(
      "02_intermediate_data/05_exposure_source_replication/GCST90018849/",
      name,
      ".csv"
    ),
    row.names = FALSE
  )
}


save_external_harm(
  mr_swedish_v_primary,
  "Swedish2026_Veillonella_primary_MR_ready"
)

save_external_harm(
  mr_swedish_va_primary,
  "Swedish2026_Veillonellaceae_primary_MR_ready"
)

save_external_harm(
  mr_dmp_v_primary,
  "DMP2022_Veillonella_primary_MR_ready"
)

save_external_harm(
  mr_dmp_va_primary,
  "DMP2022_Veillonellaceae_primary_MR_ready"
)

save_external_harm(
  mr_swedish_v_p1e5,
  "Swedish2026_Veillonella_p1e5_MR_ready"
)

save_external_harm(
  mr_swedish_va_p1e5,
  "Swedish2026_Veillonellaceae_p1e5_MR_ready"
)

save_external_harm(
  mr_dmp_v_p1e5,
  "DMP2022_Veillonella_p1e5_MR_ready"
)

save_external_harm(
  mr_dmp_va_p1e5,
  "DMP2022_Veillonellaceae_p1e5_MR_ready"
)
# ============================================================
# 42. MR analysis helper
# ============================================================

external_mr_methods <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_weighted_mode",
  "mr_simple_mode"
)


run_external_mr <- function(dat, analysis_label, role) {
  
  if (nrow(dat) < 2) {
    
    warning(
      paste(
        analysis_label,
        "has fewer than 2 SNPs; multi-SNP MR not performed."
      )
    )
    
    return(NULL)
  }
  
  res <- TwoSampleMR::mr(
    dat,
    method_list = external_mr_methods
  )
  
  res$analysis <- analysis_label
  res$analysis_role <- role
  
  # Gastric cancer is binary
  res$OR <- exp(res$b)
  res$OR_lower95 <- exp(
    res$b - 1.96 * res$se
  )
  res$OR_upper95 <- exp(
    res$b + 1.96 * res$se
  )
  
  res
}
# ============================================================
# 43. Primary external exposure-source MR
# ============================================================

mrres_swedish_v_primary <- run_external_mr(
  mr_swedish_v_primary,
  "Swedish2026_Veillonella_primary_p5e6",
  "Primary exploratory exposure-source replication"
)

mrres_swedish_va_primary <- run_external_mr(
  mr_swedish_va_primary,
  "Swedish2026_Veillonellaceae_primary_p5e6",
  "Primary exploratory exposure-source replication"
)

mrres_dmp_v_primary <- run_external_mr(
  mr_dmp_v_primary,
  "DMP2022_Veillonella_primary_p5e6",
  "Secondary exploratory exposure-source replication"
)

mrres_dmp_va_primary <- run_external_mr(
  mr_dmp_va_primary,
  "DMP2022_Veillonellaceae_primary_p5e6",
  "Secondary exploratory exposure-source replication"
)
# ============================================================
# 44. Threshold-sensitivity external MR
# ============================================================

mrres_swedish_v_p1e5 <- run_external_mr(
  mr_swedish_v_p1e5,
  "Swedish2026_Veillonella_sensitivity_p1e5",
  "Threshold sensitivity"
)

mrres_swedish_va_p1e5 <- run_external_mr(
  mr_swedish_va_p1e5,
  "Swedish2026_Veillonellaceae_sensitivity_p1e5",
  "Threshold sensitivity"
)

mrres_dmp_v_p1e5 <- run_external_mr(
  mr_dmp_v_p1e5,
  "DMP2022_Veillonella_sensitivity_p1e5",
  "Threshold sensitivity"
)

mrres_dmp_va_p1e5 <- run_external_mr(
  mr_dmp_va_p1e5,
  "DMP2022_Veillonellaceae_sensitivity_p1e5",
  "Threshold sensitivity"
)
# ============================================================
# 45. Combine all external MR results
# ============================================================

external_mr_results <- dplyr::bind_rows(
  
  mrres_swedish_v_primary,
  mrres_swedish_va_primary,
  mrres_dmp_v_primary,
  mrres_dmp_va_primary,
  
  mrres_swedish_v_p1e5,
  mrres_swedish_va_p1e5,
  mrres_dmp_v_p1e5,
  mrres_dmp_va_p1e5
)


dir.create(
  "05_results/07_exposure_source_replication",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  external_mr_results,
  "05_results/07_exposure_source_replication/all_external_MR_results.csv",
  row.names = FALSE
)
  row.names = FALSE
  # ============================================================
  # 47. All external exposure-source IVW results
  #     Primary + threshold sensitivity
  # ============================================================
  
  external_ivw_all <- external_mr_results %>%
    
    dplyr::filter(
      method == "Inverse variance weighted"
    ) %>%
    
    dplyr::mutate(
      
      dataset = dplyr::case_when(
        grepl("^Swedish2026", analysis) ~ "Swedish2026",
        grepl("^DMP2022", analysis) ~ "DMP2022",
        TRUE ~ NA_character_
      ),
      
      taxon = dplyr::case_when(
        grepl("Veillonellaceae", analysis) ~ "Veillonellaceae",
        grepl("Veillonella", analysis) ~ "Veillonella",
        TRUE ~ NA_character_
      ),
      
      threshold = dplyr::case_when(
        grepl("primary_p5e6", analysis) ~ "P <= 5e-6",
        grepl("sensitivity_p1e5", analysis) ~ "P <= 1e-5",
        TRUE ~ NA_character_
      ),
      
      multiplicity_threshold = 0.0125,
      
      significant_0.0125 = pval < multiplicity_threshold
      
    ) %>%
    
    dplyr::select(
      dataset,
      taxon,
      threshold,
      nsnp,
      b,
      se,
      OR,
      OR_lower95,
      OR_upper95,
      pval,
      significant_0.0125
    ) %>%
    
    dplyr::arrange(
      dataset,
      taxon,
      threshold
    )
  
  
  cat("\n========================================\n")
  cat("ALL EXTERNAL EXPOSURE-SOURCE IVW RESULTS\n")
  cat("========================================\n")
  
  print(
    external_ivw_all,
    row.names = FALSE
  )
#保存
  write.csv(
    external_ivw_all,
    "05_results/07_exposure_source_replication/all_external_IVW_results.csv",
    row.names = FALSE
  )
  # ============================================================
  # 48. Primary vs threshold-sensitivity comparison
  # ============================================================
  
  external_ivw_comparison <- external_ivw_all %>%
    
    dplyr::mutate(
      
      OR_95CI = sprintf(
        "%.3f (%.3f-%.3f)",
        OR,
        OR_lower95,
        OR_upper95
      ),
      
      p_display = format.pval(
        pval,
        digits = 3,
        eps = 0.001
      )
    ) %>%
    
    dplyr::select(
      dataset,
      taxon,
      threshold,
      nsnp,
      OR_95CI,
      p_display,
      significant_0.0125
    )
  
  
  cat("\n========================================\n")
  cat("PRIMARY VS THRESHOLD SENSITIVITY\n")
  cat("========================================\n")
  
  print(
    external_ivw_comparison,
    row.names = FALSE
  )
  
  
  write.csv(
    external_ivw_comparison,
    "05_results/07_exposure_source_replication/primary_vs_threshold_IVW_comparison.csv",
    row.names = FALSE
  )
  # ============================================================
  # 49. Heterogeneity and Egger-intercept helper
  # ============================================================
  
  run_external_qc_basic <- function(dat, analysis_label) {
    
    # -------------------------
    # Heterogeneity
    # -------------------------
    
    heterogeneity <- tryCatch(
      TwoSampleMR::mr_heterogeneity(dat),
      error = function(e) NULL
    )
    
    if (!is.null(heterogeneity)) {
      heterogeneity$analysis <- analysis_label
    }
    
    
    # -------------------------
    # Egger intercept
    # -------------------------
    
    egger <- tryCatch(
      TwoSampleMR::mr_pleiotropy_test(dat),
      error = function(e) NULL
    )
    
    if (!is.null(egger)) {
      egger$analysis <- analysis_label
    }
    
    
    list(
      heterogeneity = heterogeneity,
      egger = egger
    )
  }
  # ============================================================
  # 50. Run basic QC for all 8 analyses
  # ============================================================
  
  qc_swedish_v_primary <- run_external_qc_basic(
    mr_swedish_v_primary,
    "Swedish2026_Veillonella_primary_p5e6"
  )
  
  qc_swedish_va_primary <- run_external_qc_basic(
    mr_swedish_va_primary,
    "Swedish2026_Veillonellaceae_primary_p5e6"
  )
  
  qc_dmp_v_primary <- run_external_qc_basic(
    mr_dmp_v_primary,
    "DMP2022_Veillonella_primary_p5e6"
  )
  
  qc_dmp_va_primary <- run_external_qc_basic(
    mr_dmp_va_primary,
    "DMP2022_Veillonellaceae_primary_p5e6"
  )
  
  
  qc_swedish_v_p1e5 <- run_external_qc_basic(
    mr_swedish_v_p1e5,
    "Swedish2026_Veillonella_sensitivity_p1e5"
  )
  
  qc_swedish_va_p1e5 <- run_external_qc_basic(
    mr_swedish_va_p1e5,
    "Swedish2026_Veillonellaceae_sensitivity_p1e5"
  )
  
  qc_dmp_v_p1e5 <- run_external_qc_basic(
    mr_dmp_v_p1e5,
    "DMP2022_Veillonella_sensitivity_p1e5"
  )
  
  qc_dmp_va_p1e5 <- run_external_qc_basic(
    mr_dmp_va_p1e5,
    "DMP2022_Veillonellaceae_sensitivity_p1e5"
  )
  # ============================================================
  # 51. Combine heterogeneity results
  # ============================================================
  
  external_heterogeneity <- dplyr::bind_rows(
    
    qc_swedish_v_primary$heterogeneity,
    qc_swedish_va_primary$heterogeneity,
    qc_dmp_v_primary$heterogeneity,
    qc_dmp_va_primary$heterogeneity,
    
    qc_swedish_v_p1e5$heterogeneity,
    qc_swedish_va_p1e5$heterogeneity,
    qc_dmp_v_p1e5$heterogeneity,
    qc_dmp_va_p1e5$heterogeneity
  )
  
  
  cat("\n========================================\n")
  cat("EXTERNAL MR HETEROGENEITY\n")
  cat("========================================\n")
  
  print(
    external_heterogeneity,
    row.names = FALSE
  )
  # ============================================================
  # 52. Combine Egger-intercept results
  # ============================================================
  
  external_egger <- dplyr::bind_rows(
    
    qc_swedish_v_primary$egger,
    qc_swedish_va_primary$egger,
    qc_dmp_v_primary$egger,
    qc_dmp_va_primary$egger,
    
    qc_swedish_v_p1e5$egger,
    qc_swedish_va_p1e5$egger,
    qc_dmp_v_p1e5$egger,
    qc_dmp_va_p1e5$egger
  )
  
  
  cat("\n========================================\n")
  cat("EXTERNAL MR EGGER INTERCEPT\n")
  cat("========================================\n")
  
  print(
    external_egger,
    row.names = FALSE
  )
  #保存
  write.csv(
    external_heterogeneity,
    "05_results/07_exposure_source_replication/external_heterogeneity.csv",
    row.names = FALSE
  )
  
  write.csv(
    external_egger,
    "05_results/07_exposure_source_replication/external_Egger_intercept.csv",
    row.names = FALSE
  )
  
  # ============================================================
  # 53. MR-PRESSO helper
  # ============================================================
  
  run_external_presso <- function(dat, analysis_label) {
    
    if (nrow(dat) < 4) {
      
      cat(
        "\nMR-PRESSO skipped:",
        analysis_label,
        "- fewer than 4 SNPs\n"
      )
      
      return(
        list(
          analysis = analysis_label,
          nsnp = nrow(dat),
          feasible = FALSE,
          result = NULL
        )
      )
    }
    
    cat(
      "\nRunning MR-PRESSO:",
      analysis_label,
      "with",
      nrow(dat),
      "SNPs\n"
    )
    
    res <- tryCatch(
      
      MRPRESSO::mr_presso(
        
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        
        data = dat,
        
        NbDistribution = 10000,
        SignifThreshold = 0.05
      ),
      
      error = function(e) {
        
        message(
          "MR-PRESSO error for ",
          analysis_label,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    list(
      analysis = analysis_label,
      nsnp = nrow(dat),
      feasible = !is.null(res),
      result = res
    )
  }
  # ============================================================
  # 54. Run MR-PRESSO for all external analyses
  # ============================================================
  
  presso_swedish_v_primary <- run_external_presso(
    mr_swedish_v_primary,
    "Swedish2026_Veillonella_primary_p5e6"
  )
  
  presso_swedish_va_primary <- run_external_presso(
    mr_swedish_va_primary,
    "Swedish2026_Veillonellaceae_primary_p5e6"
  )
  
  presso_dmp_v_primary <- run_external_presso(
    mr_dmp_v_primary,
    "DMP2022_Veillonella_primary_p5e6"
  )
  
  presso_dmp_va_primary <- run_external_presso(
    mr_dmp_va_primary,
    "DMP2022_Veillonellaceae_primary_p5e6"
  )
  
  
  presso_swedish_v_p1e5 <- run_external_presso(
    mr_swedish_v_p1e5,
    "Swedish2026_Veillonella_sensitivity_p1e5"
  )
  
  presso_swedish_va_p1e5 <- run_external_presso(
    mr_swedish_va_p1e5,
    "Swedish2026_Veillonellaceae_sensitivity_p1e5"
  )
  
  presso_dmp_v_p1e5 <- run_external_presso(
    mr_dmp_v_p1e5,
    "DMP2022_Veillonella_sensitivity_p1e5"
  )
  
  presso_dmp_va_p1e5 <- run_external_presso(
    mr_dmp_va_p1e5,
    "DMP2022_Veillonellaceae_sensitivity_p1e5"
  )
  # ============================================================
  # 55. Extract MR-PRESSO global-test P values
  # ============================================================
  
  extract_presso_global <- function(presso_object) {
    
    if (
      is.null(presso_object$result) ||
      !presso_object$feasible
    ) {
      
      return(
        data.frame(
          analysis = presso_object$analysis,
          nsnp = presso_object$nsnp,
          feasible = FALSE,
          global_p = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    res <- presso_object$result
    
    global_p <- tryCatch(
      
      res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
      
      error = function(e) NA
    )
    
    global_p <- as.numeric(
      sub(
        "^<",
        "",
        as.character(global_p)
      )
    )
    
    data.frame(
      analysis = presso_object$analysis,
      nsnp = presso_object$nsnp,
      feasible = TRUE,
      global_p = global_p,
      stringsAsFactors = FALSE
    )
  }
  #汇总
  external_presso_summary <- dplyr::bind_rows(
    
    extract_presso_global(presso_swedish_v_primary),
    extract_presso_global(presso_swedish_va_primary),
    extract_presso_global(presso_dmp_v_primary),
    extract_presso_global(presso_dmp_va_primary),
    
    extract_presso_global(presso_swedish_v_p1e5),
    extract_presso_global(presso_swedish_va_p1e5),
    extract_presso_global(presso_dmp_v_p1e5),
    extract_presso_global(presso_dmp_va_p1e5)
  )
  
  
  cat("\n========================================\n")
  cat("EXTERNAL MR-PRESSO GLOBAL TEST\n")
  cat("========================================\n")
  
  print(
    external_presso_summary,
    row.names = FALSE
  )
  # ============================================================
  # 56. Save MR-PRESSO objects
  # ============================================================
  
  saveRDS(
    list(
      Swedish_V_primary = presso_swedish_v_primary,
      Swedish_Va_primary = presso_swedish_va_primary,
      DMP_V_primary = presso_dmp_v_primary,
      DMP_Va_primary = presso_dmp_va_primary,
      
      Swedish_V_p1e5 = presso_swedish_v_p1e5,
      Swedish_Va_p1e5 = presso_swedish_va_p1e5,
      DMP_V_p1e5 = presso_dmp_v_p1e5,
      DMP_Va_p1e5 = presso_dmp_va_p1e5
    ),
    
    "05_results/07_exposure_source_replication/external_MRPRESSO_all.rds"
  )
  
  
  write.csv(
    external_presso_summary,
    "05_results/07_exposure_source_replication/external_MRPRESSO_summary.csv",
    row.names = FALSE
  )
  # ============================================================
  # 57. Leave-one-out and single-SNP helper
  # ============================================================
  
  run_external_influence <- function(dat, analysis_label) {
    
    loo <- tryCatch(
      TwoSampleMR::mr_leaveoneout(dat),
      error = function(e) NULL
    )
    
    single <- tryCatch(
      TwoSampleMR::mr_singlesnp(dat),
      error = function(e) NULL
    )
    
    if (!is.null(loo)) {
      loo$analysis <- analysis_label
    }
    
    if (!is.null(single)) {
      single$analysis <- analysis_label
    }
    
    list(
      leaveoneout = loo,
      singlesnp = single
    )
  }
  # ============================================================
  # 58. Run influence diagnostics
  # ============================================================
  
  inf_swedish_v_primary <- run_external_influence(
    mr_swedish_v_primary,
    "Swedish2026_Veillonella_primary_p5e6"
  )
  
  inf_swedish_va_primary <- run_external_influence(
    mr_swedish_va_primary,
    "Swedish2026_Veillonellaceae_primary_p5e6"
  )
  
  inf_dmp_v_primary <- run_external_influence(
    mr_dmp_v_primary,
    "DMP2022_Veillonella_primary_p5e6"
  )
  
  inf_dmp_va_primary <- run_external_influence(
    mr_dmp_va_primary,
    "DMP2022_Veillonellaceae_primary_p5e6"
  )
  
  
  inf_swedish_v_p1e5 <- run_external_influence(
    mr_swedish_v_p1e5,
    "Swedish2026_Veillonella_sensitivity_p1e5"
  )
  
  inf_swedish_va_p1e5 <- run_external_influence(
    mr_swedish_va_p1e5,
    "Swedish2026_Veillonellaceae_sensitivity_p1e5"
  )
  
  inf_dmp_v_p1e5 <- run_external_influence(
    mr_dmp_v_p1e5,
    "DMP2022_Veillonella_sensitivity_p1e5"
  )
  
  inf_dmp_va_p1e5 <- run_external_influence(
    mr_dmp_va_p1e5,
    "DMP2022_Veillonellaceae_sensitivity_p1e5"
  )
  # ============================================================
  # 59. Combine/save LOO and single-SNP
  # ============================================================
  
  external_leaveoneout <- dplyr::bind_rows(
    inf_swedish_v_primary$leaveoneout,
    inf_swedish_va_primary$leaveoneout,
    inf_dmp_v_primary$leaveoneout,
    inf_dmp_va_primary$leaveoneout,
    
    inf_swedish_v_p1e5$leaveoneout,
    inf_swedish_va_p1e5$leaveoneout,
    inf_dmp_v_p1e5$leaveoneout,
    inf_dmp_va_p1e5$leaveoneout
  )
  
  
  external_singlesnp <- dplyr::bind_rows(
    inf_swedish_v_primary$singlesnp,
    inf_swedish_va_primary$singlesnp,
    inf_dmp_v_primary$singlesnp,
    inf_dmp_va_primary$singlesnp,
    
    inf_swedish_v_p1e5$singlesnp,
    inf_swedish_va_p1e5$singlesnp,
    inf_dmp_v_p1e5$singlesnp,
    inf_dmp_va_p1e5$singlesnp
  )
  
  
  write.csv(
    external_leaveoneout,
    "05_results/07_exposure_source_replication/external_leaveoneout.csv",
    row.names = FALSE
  )
  
  write.csv(
    external_singlesnp,
    "05_results/07_exposure_source_replication/external_singleSNP.csv",
    row.names = FALSE
  )
  # ============================================================
  # 60. Inspect key DMP Veillonellaceae sensitivity analysis
  # ============================================================
  
  key_loo <- external_leaveoneout[
    external_leaveoneout$analysis ==
      "DMP2022_Veillonellaceae_sensitivity_p1e5",
  ]
  
  key_single <- external_singlesnp[
    external_singlesnp$analysis ==
      "DMP2022_Veillonellaceae_sensitivity_p1e5",
  ]
  
  cat("\n========================================\n")
  cat("KEY LOO: DMP VEILLONELLACEAE P1E-5\n")
  cat("========================================\n")
  
  print(
    key_loo,
    row.names = FALSE
  )
  
  cat("\n========================================\n")
  cat("KEY SINGLE-SNP: DMP VEILLONELLACEAE P1E-5\n")
  cat("========================================\n")
  
  print(
    key_single,
    row.names = FALSE
  )
  # ============================================================
  # 61. Leave-one-out summary for all external analyses
  # ============================================================
  
  external_loo_summary <- external_leaveoneout %>%
    
    dplyr::filter(
      !grepl("^All", SNP)
    ) %>%
    
    dplyr::group_by(analysis) %>%
    
    dplyr::summarise(
      n_leave_one_out = dplyr::n(),
      all_beta_positive = all(b > 0, na.rm = TRUE),
      all_beta_negative = all(b < 0, na.rm = TRUE),
      min_loo_beta = min(b, na.rm = TRUE),
      max_loo_beta = max(b, na.rm = TRUE),
      min_loo_p = min(p, na.rm = TRUE),
      max_loo_p = max(p, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  print(
    external_loo_summary,
    row.names = FALSE
  )
  # ============================================================
  # 62. Build final external-replication QC summary
  # ============================================================
  
  external_ivw_qc_base <- external_mr_results %>%
    
    dplyr::filter(
      method == "Inverse variance weighted"
    ) %>%
    
    dplyr::select(
      analysis,
      nsnp,
      b,
      se,
      OR,
      OR_lower95,
      OR_upper95,
      pval
    )
  
  
  external_het_ivw <- external_heterogeneity %>%
    
    dplyr::filter(
      method == "Inverse variance weighted"
    ) %>%
    
    dplyr::select(
      analysis,
      Q,
      Q_df,
      Q_pval
    ) %>%
    
    dplyr::rename(
      IVW_Q = Q,
      IVW_Q_df = Q_df,
      IVW_Q_p = Q_pval
    )
  
  
  external_egger_qc <- external_egger %>%
    
    dplyr::select(
      analysis,
      egger_intercept,
      se,
      pval
    ) %>%
    
    dplyr::rename(
      Egger_intercept = egger_intercept,
      Egger_intercept_SE = se,
      Egger_intercept_p = pval
    )
  
  
  external_presso_qc <- external_presso_summary %>%
    
    dplyr::select(
      analysis,
      feasible,
      global_p
    ) %>%
    
    dplyr::rename(
      MRPRESSO_feasible = feasible,
      MRPRESSO_global_p = global_p
    )
  
  
  external_final_qc <- external_ivw_qc_base %>%
    
    dplyr::left_join(
      external_het_ivw,
      by = "analysis"
    ) %>%
    
    dplyr::left_join(
      external_egger_qc,
      by = "analysis"
    ) %>%
    
    dplyr::left_join(
      external_presso_qc,
      by = "analysis"
    ) %>%
    
    dplyr::left_join(
      external_loo_summary,
      by = "analysis"
    ) %>%
    
    dplyr::mutate(
      significant_0.0125 = pval < 0.0125,
      
      heterogeneity_evidence =
        !is.na(IVW_Q_p) & IVW_Q_p < 0.05,
      
      directional_pleiotropy_evidence =
        !is.na(Egger_intercept_p) &
        Egger_intercept_p < 0.05,
      
      presso_global_evidence =
        !is.na(MRPRESSO_global_p) &
        MRPRESSO_global_p < 0.05
    )
  #打印最重要的列
  cat("\n========================================\n")
  cat("FINAL EXTERNAL REPLICATION QC SUMMARY\n")
  cat("========================================\n")
  
  print(
    external_final_qc[
      ,
      c(
        "analysis",
        "nsnp",
        "OR",
        "OR_lower95",
        "OR_upper95",
        "pval",
        "IVW_Q_p",
        "Egger_intercept_p",
        "MRPRESSO_global_p",
        "all_beta_positive",
        "all_beta_negative"
      )
    ],
    row.names = FALSE
  )
  #保存
  write.csv(
    external_final_qc,
    "05_results/07_exposure_source_replication/external_replication_final_QC_summary.csv",
    row.names = FALSE
  )
  #检查
  cat(
    "\nHeterogeneity evidence:",
    sum(external_final_qc$heterogeneity_evidence, na.rm = TRUE),
    "/ 8\n"
  )
  
  cat(
    "Directional pleiotropy evidence:",
    sum(external_final_qc$directional_pleiotropy_evidence, na.rm = TRUE),
    "/ 8\n"
  )
  
  cat(
    "MR-PRESSO global evidence:",
    sum(external_final_qc$presso_global_evidence, na.rm = TRUE),
    "/ 8",
    "(one analysis may be not feasible)\n"
  )
  # ============================================================
  # 63. Load frozen MiBioGen core IVW results
  # ============================================================
  
  core_ivw <- read.csv(
    "05_results/06_summary/IVW_summary.csv",
    stringsAsFactors = FALSE
  )
  
  cat("\nCore IVW columns:\n")
  print(names(core_ivw))
  # ============================================================
  # 64. Extract MiBioGen -> GCST90018849
  #     primary + threshold sensitivity
  # ============================================================
  
  mibiogen_exposure_source <- core_ivw %>%
    
    dplyr::filter(
      direction == "Forward",
      outcome_label == "Gastric cancer",
      exposure_label %in% c(
        "Veillonella",
        "Veillonellaceae"
      ),
      analysis %in% c(
        "Primary GCST90018849",
        "Threshold sensitivity P<=1e-5"
      )
    ) %>%
    
    dplyr::mutate(
      
      source = "MiBioGen",
      
      threshold = dplyr::case_when(
        analysis == "Primary GCST90018849" ~ "P <= 5e-6",
        grepl("Threshold", analysis) ~ "P <= 1e-5",
        TRUE ~ NA_character_
      ),
      
      taxon = exposure_label,
      
      OR = estimate,
      OR_lower95 = ci_lower,
      OR_upper95 = ci_upper
      
    ) %>%
    
    dplyr::select(
      source,
      taxon,
      threshold,
      nsnp,
      OR,
      OR_lower95,
      OR_upper95,
      pval
    )
  
  
  cat("\n========================================\n")
  cat("MIBIOGEN EXPOSURE-SOURCE ROWS\n")
  cat("========================================\n")
  
  print(
    mibiogen_exposure_source,
    row.names = FALSE
  )
  # ============================================================
  # 65. External exposure-source rows
  # ============================================================
  
  external_exposure_source <- external_ivw_all %>%
    
    dplyr::mutate(
      source = dplyr::case_when(
        dataset == "Swedish2026" ~ "Swedish 2026",
        dataset == "DMP2022" ~ "DMP 2022",
        TRUE ~ dataset
      )
    ) %>%
    
    dplyr::select(
      source,
      taxon,
      threshold,
      nsnp,
      OR,
      OR_lower95,
      OR_upper95,
      pval
    )
  # ============================================================
  # 66. Exposure-source robustness summary
  # ============================================================
  
  exposure_source_robustness <- dplyr::bind_rows(
    mibiogen_exposure_source,
    external_exposure_source
  ) %>%
    
    dplyr::mutate(
      
      source = factor(
        source,
        levels = c(
          "MiBioGen",
          "Swedish 2026",
          "DMP 2022"
        )
      ),
      
      threshold = factor(
        threshold,
        levels = c(
          "P <= 5e-6",
          "P <= 1e-5"
        )
      ),
      
      significant_0.0125 =
        pval < 0.0125,
      
      OR_95CI = sprintf(
        "%.3f (%.3f-%.3f)",
        OR,
        OR_lower95,
        OR_upper95
      ),
      
      p_display = format.pval(
        pval,
        digits = 3,
        eps = 0.001
      )
    ) %>%
    
    dplyr::arrange(
      taxon,
      source,
      threshold
    )
  
  
  cat("\n========================================\n")
  cat("EXPOSURE-SOURCE ROBUSTNESS SUMMARY\n")
  cat("========================================\n")
  
  print(
    exposure_source_robustness[
      ,
      c(
        "source",
        "taxon",
        "threshold",
        "nsnp",
        "OR_95CI",
        "p_display",
        "significant_0.0125"
      )
    ],
    row.names = FALSE
  )
  
  
  write.csv(
    exposure_source_robustness,
    "05_results/07_exposure_source_replication/exposure_source_robustness_summary.csv",
    row.names = FALSE
  )
  # ============================================================
  # 67. Exposure-source robustness forest plot
  # ============================================================
  
  dir.create(
    "06_figures/07_exposure_source_replication",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  plot_dat <- exposure_source_robustness %>%
    
    dplyr::mutate(
      
      analysis_label = paste(
        source,
        threshold,
        sep = " | "
      ),
      
      analysis_label = factor(
        analysis_label,
        levels = rev(
          c(
            "MiBioGen | P <= 5e-6",
            "MiBioGen | P <= 1e-5",
            "Swedish 2026 | P <= 5e-6",
            "Swedish 2026 | P <= 1e-5",
            "DMP 2022 | P <= 5e-6",
            "DMP 2022 | P <= 1e-5"
          )
        )
      )
    )
  
  
  exposure_source_forest <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = OR,
      y = analysis_label
    )
  ) +
    
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = OR_lower95,
        xmax = OR_upper95
      ),
      orientation = "y",
      width = 0.15,
      linewidth = 0.7
    ) +
    
    ggplot2::geom_point(
      size = 2.8
    ) +
    
    ggplot2::facet_wrap(
      ~ taxon,
      ncol = 1,
      scales = "free_y"
    ) +
    
    ggplot2::scale_x_log10(
      breaks = c(
        0.6,
        0.8,
        1,
        1.2,
        1.5,
        2
      )
    ) +
    
    ggplot2::labs(
      x = "Odds ratio for gastric cancer (95% CI)",
      y = NULL,
      title = NULL
    ) +
    
    ggplot2::theme_bw(
      base_size = 12
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  
  exposure_source_forest
 #保存
  ggplot2::ggsave(
    "06_figures/07_exposure_source_replication/Exposure_source_robustness_forest.pdf",
    plot = exposure_source_forest,
    width = 8,
    height = 7
  )
  
  ggplot2::ggsave(
    "06_figures/07_exposure_source_replication/Exposure_source_robustness_forest.png",
    plot = exposure_source_forest,
    width = 8,
    height = 7,
    dpi = 600
  )
  #检查
  cat("\nRows in exposure_source_robustness:\n")
  print(nrow(exposure_source_robustness))
  
  cat("\n========================================\n")
  cat("EXPOSURE-SOURCE ROBUSTNESS SUMMARY\n")
  cat("========================================\n")
  
  print(
    exposure_source_robustness[
      ,
      c(
        "source",
        "taxon",
        "threshold",
        "nsnp",
        "OR_95CI",
        "p_display",
        "significant_0.0125"
      )
    ],
    row.names = FALSE
  )
  #检查文件
  cat("\nResults files:\n")
  print(
    list.files(
      "05_results/07_exposure_source_replication"
    )
  )
  
  cat("\nFigure files:\n")
  print(
    list.files(
      "06_figures/07_exposure_source_replication"
    )
  )
 
  
  # ============================================================
  # 69. Final post-registration exposure-replication audit
  # ============================================================
  
  final_audit_text <- c(
    
    "PROJECT",
    "01_MR_GastricCancer",
    "",
    "DOCUMENT",
    "Post-registration Exposure-Source Replication Final Audit v1.0",
    "",
    "DATE",
    "7 September 2026",
    "",
    "OSF REGISTRATION",
    "Secondary Data Preregistration",
    "Registration ID: Kv7tj",
    "",
    "The exploratory exposure-source replication was separately",
    "registered before generation or inspection of the corresponding",
    "Swedish 2026 / DMP -> GCST90018849 MR results.",
    "",
    "============================================================",
    "1. REGISTERED EXPOSURE SOURCES",
    "============================================================",
    "",
    "Swedish 2026:",
    "Veillonella       GCST90671339",
    "Veillonellaceae   GCST90671681",
    "",
    "Dutch Microbiome Project:",
    "Veillonella       GCST90027724",
    "Veillonellaceae   GCST90027679",
    "",
    "Fixed gastric cancer outcome:",
    "GCST90018849",
    "",
    "============================================================",
    "2. PRIMARY INSTRUMENT SETS",
    "============================================================",
    "",
    "Primary exposure instrument threshold:",
    "P <= 5e-6",
    "",
    "LD clumping:",
    "r2 < 0.001",
    "10,000 kb",
    "EUR reference",
    "",
    "Instrument strength criterion:",
    "F > 10",
    "",
    "Final clumped exposure IVs:",
    "",
    "Swedish2026 Veillonella:       7",
    "Swedish2026 Veillonellaceae:   8",
    "DMP2022 Veillonella:           5",
    "DMP2022 Veillonellaceae:       7",
    "",
    "After conservative action=3 harmonisation:",
    "",
    "Swedish2026 Veillonella:       4",
    "Swedish2026 Veillonellaceae:   3",
    "DMP2022 Veillonella:           4",
    "DMP2022 Veillonellaceae:       6",
    "",
    "============================================================",
    "3. THRESHOLD-SENSITIVITY INSTRUMENT SETS",
    "============================================================",
    "",
    "Prespecified secondary sensitivity threshold:",
    "P <= 1e-5",
    "",
    "Final clumped exposure IVs:",
    "",
    "Swedish2026 Veillonella:       15",
    "Swedish2026 Veillonellaceae:   18",
    "DMP2022 Veillonella:            8",
    "DMP2022 Veillonellaceae:       13",
    "",
    "After harmonisation:",
    "",
    "Swedish2026 Veillonella:       10",
    "Swedish2026 Veillonellaceae:   10",
    "DMP2022 Veillonella:            6",
    "DMP2022 Veillonellaceae:       12",
    "",
    "============================================================",
    "4. OUTCOME EXTRACTION",
    "============================================================",
    "",
    "Unique external IV SNPs requested from GCST90018849:",
    "42",
    "",
    "Found:",
    "42",
    "",
    "Missing:",
    "0",
    "",
    "Missing beta:",
    "0",
    "",
    "Missing standard error:",
    "0",
    "",
    "Missing P value:",
    "0",
    "",
    "Duplicate SNPs:",
    "0",
    "",
    "No proxy variant was required.",
    "",
    "============================================================",
    "5. PRIMARY EXTERNAL IVW RESULTS",
    "============================================================",
    "",
    "Swedish2026 Veillonella:",
    "OR 0.921",
    "95% CI 0.653-1.300",
    "P = 0.639",
    "",
    "Swedish2026 Veillonellaceae:",
    "OR 0.878",
    "95% CI 0.601-1.283",
    "P = 0.503",
    "",
    "DMP2022 Veillonella:",
    "OR 1.018",
    "95% CI 0.901-1.150",
    "P = 0.775",
    "",
    "DMP2022 Veillonellaceae:",
    "OR 1.145",
    "95% CI 0.934-1.403",
    "P = 0.192",
    "",
    "None of the four primary external comparisons met the",
    "multiplicity-aware reference threshold P < 0.0125.",
    "",
    "============================================================",
    "6. THRESHOLD-SENSITIVITY RESULTS",
    "============================================================",
    "",
    "Swedish2026 Veillonella:",
    "OR 0.915",
    "95% CI 0.684-1.225",
    "P = 0.553",
    "",
    "Swedish2026 Veillonellaceae:",
    "OR 0.925",
    "95% CI 0.621-1.377",
    "P = 0.700",
    "",
    "DMP2022 Veillonella:",
    "OR 1.040",
    "95% CI 0.933-1.160",
    "P = 0.479",
    "",
    "DMP2022 Veillonellaceae:",
    "OR 1.167",
    "95% CI 1.045-1.303",
    "P = 0.00616",
    "",
    "The DMP2022 Veillonellaceae P <= 1e-5 analysis met the",
    "multiplicity-aware reference threshold.",
    "",
    "This remains a prespecified secondary threshold-sensitivity",
    "finding and does not replace the nonsignificant primary",
    "P <= 5e-6 analysis.",
    "",
    "============================================================",
    "7. QC SUMMARY",
    "============================================================",
    "",
    "Heterogeneity evidence:",
    "0 / 8 analyses",
    "",
    "Directional horizontal pleiotropy evidence by MR-Egger intercept:",
    "0 / 8 analyses",
    "",
    "MR-PRESSO global evidence:",
    "0 / 8 analyses",
    "",
    "MR-PRESSO was not methodologically feasible for the",
    "3-SNP Swedish2026 Veillonellaceae primary analysis.",
    "",
    "For DMP2022 Veillonellaceae P <= 1e-5:",
    "",
    "IVW heterogeneity P = 0.718",
    "MR-Egger intercept P = 0.871",
    "MR-PRESSO global P = 0.7752",
    "No MR-PRESSO outlier was identified.",
    "",
    "Leave-one-out analysis retained a positive overall beta",
    "after exclusion of every individual genetic instrument.",
    "",
    "============================================================",
    "8. EXPOSURE-SOURCE ROBUSTNESS",
    "============================================================",
    "",
    "Veillonella external exposure-source estimates were generally",
    "compatible with the null.",
    "",
    "For Veillonellaceae, MiBioGen and DMP estimates were",
    "directionally positive, whereas Swedish 2026 estimates were",
    "negative at both instrument thresholds.",
    "",
    "The direction of the Veillonellaceae MR estimate was therefore",
    "not consistently robust to the microbiome exposure GWAS source.",
    "",
    "No formal meta-analysis across MiBioGen, Swedish 2026, and DMP",
    "was performed because microbial phenotype scales and preprocessing",
    "procedures were not assumed to be directly commensurate across",
    "GWAS sources.",
    "",
    "============================================================",
    "9. INTERPRETATION LOCK",
    "============================================================",
    "",
    "The post-registration exploratory analyses do not establish",
    "a uniformly replicated causal association.",
    "",
    "The significant DMP threshold-sensitivity result will be",
    "reported transparently as a secondary sensitivity finding",
    "and will not be used to retrospectively redefine the primary",
    "analysis.",
    "",
    "Absence of statistically significant heterogeneity, MR-Egger",
    "intercept, or MR-PRESSO global tests will be reported as",
    "'no statistical evidence of' the corresponding violation and",
    "will not be interpreted as proof of absence.",
    "",
    "============================================================",
    "10. ANALYTICAL STATUS",
    "============================================================",
    "",
    "The registered post-registration exploratory exposure-source",
    "replication analysis is complete.",
    "",
    "No analytical rule was changed on the basis of observed MR",
    "effect estimates, confidence intervals, or P values.",
    "",
    "The analytical outputs documented in this audit are frozen",
    "for manuscript preparation."
  )
  
  
  writeLines(
    final_audit_text,
    "00_admin/13_postregistration_exposure_replication_final_audit_v1.0.txt"
  )
  
  
  cat("\nSaved final audit:\n")
  cat(
    "00_admin/13_postregistration_exposure_replication_final_audit_v1.0.txt\n"
  )
  
  stopifnot(
    file.exists(
      "00_admin/13_postregistration_exposure_replication_final_audit_v1.0.txt"
    )
  )
  
  # ============================================================
  # 70. Snapshot configuration
  # ============================================================
  
  snapshot_root <-
    "09_snapshot/2026-09-07_postregistration_exposure_replication_v1.0"
  
  
  # Protect against accidental overwrite
  if (dir.exists(snapshot_root)) {
    
    stop(
      paste0(
        "Snapshot already exists:\n",
        snapshot_root,
        "\nDelete/rename it manually only if you intentionally want to rebuild it."
      )
    )
  }
  
  
  dir.create(
    snapshot_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ============================================================
  # 71. Recursive copy helper
  # ============================================================
  
  copy_tree <- function(source, destination) {
    
    if (!file.exists(source) && !dir.exists(source)) {
      stop(
        paste(
          "Source does not exist:",
          source
        )
      )
    }
    
    
    # ----------------------------------------------------------
    # Single file
    # ----------------------------------------------------------
    
    if (!dir.exists(source)) {
      
      dir.create(
        dirname(destination),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      ok <- file.copy(
        source,
        destination,
        overwrite = FALSE,
        copy.mode = TRUE,
        copy.date = TRUE
      )
      
      if (!ok) {
        stop(
          paste(
            "Failed to copy:",
            source
          )
        )
      }
      
      return(invisible(TRUE))
    }
    
    
    # ----------------------------------------------------------
    # Directory tree
    # ----------------------------------------------------------
    
    files <- list.files(
      source,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE
    )
    
    dir.create(
      destination,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    if (length(files) == 0) {
      return(invisible(TRUE))
    }
    
    
    source_norm <- normalizePath(
      source,
      winslash = "/",
      mustWork = TRUE
    )
    
    files_norm <- normalizePath(
      files,
      winslash = "/",
      mustWork = TRUE
    )
    
    
    rel <- substring(
      files_norm,
      nchar(source_norm) + 2
    )
    
    
    for (i in seq_along(files)) {
      
      dest_i <- file.path(
        destination,
        rel[i]
      )
      
      
      if (dir.exists(files[i])) {
        
        dir.create(
          dest_i,
          recursive = TRUE,
          showWarnings = FALSE
        )
        
      } else {
        
        dir.create(
          dirname(dest_i),
          recursive = TRUE,
          showWarnings = FALSE
        )
        
        ok <- file.copy(
          files[i],
          dest_i,
          overwrite = FALSE,
          copy.mode = TRUE,
          copy.date = TRUE
        )
        
        if (!ok) {
          stop(
            paste(
              "Failed to copy:",
              files[i]
            )
          )
        }
      }
    }
    
    
    invisible(TRUE)
  }
  # ============================================================
  # 72. Copy administrative audit records
  # ============================================================
  
  admin_files <- c(
    
    "00_admin/10_postregistration_exposure_source_audit_v1.0.txt",
    
    "00_admin/11_postregistration_exposure_replication_plan_v1.0.txt",
    
    "00_admin/12_postregistration_harmonisation_clarification_v1.0.txt",
    
    "00_admin/13_postregistration_exposure_replication_final_audit_v1.0.txt"
  )
  
  
  stopifnot(
    all(file.exists(admin_files))
  )
  
  
  dir.create(
    file.path(snapshot_root, "00_admin"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (x in admin_files) {
    
    file.copy(
      x,
      file.path(
        snapshot_root,
        "00_admin",
        basename(x)
      ),
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
  }
  # ============================================================
  # 73. Copy exploratory intermediate datasets
  # ============================================================
  
  copy_tree(
    "02_intermediate_data/05_exposure_source_replication",
    file.path(
      snapshot_root,
      "02_intermediate_data/05_exposure_source_replication"
    )
  )
  # ============================================================
  # 74. Copy analysis script
  # ============================================================
  
  dir.create(
    file.path(snapshot_root, "04_R_scripts"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  file.copy(
    "04_R_scripts/10_exploratory_exposure_replication.R",
    file.path(
      snapshot_root,
      "04_R_scripts/10_exploratory_exposure_replication.R"
    ),
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  # ============================================================
  # 75. Copy final result outputs
  # ============================================================
  
  copy_tree(
    "05_results/07_exposure_source_replication",
    file.path(
      snapshot_root,
      "05_results/07_exposure_source_replication"
    )
  )
  # ============================================================
  # 76. Copy final figures
  # ============================================================
  
  copy_tree(
    "06_figures/07_exposure_source_replication",
    file.path(
      snapshot_root,
      "06_figures/07_exposure_source_replication"
    )
  )
  # ============================================================
  # 77. Snapshot README
  # ============================================================
  
  snapshot_readme <- c(
    
    "01_MR_GastricCancer",
    "Post-registration Exposure-Source Replication Snapshot",
    "",
    "Snapshot version:",
    "2026-09-07_postregistration_exposure_replication_v1.0",
    "",
    "Created:",
    "7 September 2026",
    "",
    "OSF Secondary Data Preregistration:",
    "Kv7tj",
    "",
    "PURPOSE",
    "This snapshot freezes the completed registered post-registration",
    "exploratory exposure-source replication analyses using Swedish 2026",
    "and Dutch Microbiome Project microbiome GWAS datasets with",
    "GCST90018849 as the fixed gastric cancer outcome.",
    "",
    "CONTENTS",
    "",
    "00_admin/",
    "Administrative audit, analysis plan, harmonisation clarification,",
    "and final analysis audit.",
    "",
    "02_intermediate_data/05_exposure_source_replication/",
    "Instrument sets, outcome extraction, harmonisation and MR-ready",
    "intermediate analytical datasets.",
    "",
    "04_R_scripts/",
    "Final exploratory exposure-source replication R script.",
    "",
    "05_results/07_exposure_source_replication/",
    "MR estimates and sensitivity/QC results.",
    "",
    "06_figures/07_exposure_source_replication/",
    "Final frozen exposure-source robustness figures.",
    "",
    "RAW DATA",
    "Large original GWAS summary-statistic files are intentionally not",
    "duplicated in this snapshot. Their accession numbers and source",
    "metadata are documented in the administrative audit and",
    "preregistration records.",
    "",
    "ANALYTICAL STATUS",
    "The registered exploratory analysis is complete.",
    "",
    "The snapshot must not be modified in place.",
    "",
    "Future manuscript-specific reformatting or graphical modifications",
    "should be generated outside this frozen analytical snapshot.",
    "",
    "Integrity is documented in snapshot_file_manifest.csv and",
    "snapshot_MD5_manifest.csv."
  )
  
  
  writeLines(
    snapshot_readme,
    file.path(
      snapshot_root,
      "README_snapshot.txt"
    )
  )
  # ============================================================
  # 78. Snapshot file manifest
  # ============================================================
  
  snapshot_files <- list.files(
    snapshot_root,
    recursive = TRUE,
    full.names = TRUE
  )
  
  snapshot_files <- snapshot_files[
    !dir.exists(snapshot_files)
  ]
  
  
  snapshot_info <- file.info(
    snapshot_files
  )
  
  
  snapshot_root_norm <- normalizePath(
    snapshot_root,
    winslash = "/",
    mustWork = TRUE
  )
  
  snapshot_files_norm <- normalizePath(
    snapshot_files,
    winslash = "/",
    mustWork = TRUE
  )
  
  
  snapshot_relative_paths <- substring(
    snapshot_files_norm,
    nchar(snapshot_root_norm) + 2
  )
  
  
  snapshot_file_manifest <- data.frame(
    
    relative_path = snapshot_relative_paths,
    
    size_bytes = snapshot_info$size,
    
    modified_time =
      as.character(
        snapshot_info$mtime
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    snapshot_file_manifest,
    file.path(
      snapshot_root,
      "snapshot_file_manifest.csv"
    ),
    row.names = FALSE
  )
  # ============================================================
  # 79. MD5 integrity manifest
  # ============================================================
  
  # Re-list because the file manifest itself now exists
  snapshot_files_md5 <- list.files(
    snapshot_root,
    recursive = TRUE,
    full.names = TRUE
  )
  
  snapshot_files_md5 <- snapshot_files_md5[
    !dir.exists(snapshot_files_md5)
  ]
  
  
  md5_values <- tools::md5sum(
    snapshot_files_md5
  )
  
  
  snapshot_files_md5_norm <- normalizePath(
    snapshot_files_md5,
    winslash = "/",
    mustWork = TRUE
  )
  
  
  md5_relative_paths <- substring(
    snapshot_files_md5_norm,
    nchar(snapshot_root_norm) + 2
  )
  
  
  snapshot_MD5_manifest <- data.frame(
    
    relative_path = md5_relative_paths,
    
    md5 = unname(md5_values),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    snapshot_MD5_manifest,
    file.path(
      snapshot_root,
      "snapshot_MD5_manifest.csv"
    ),
    row.names = FALSE
  )
  # ============================================================
  # 80. Final snapshot validation
  # ============================================================
  
  cat("\n========================================\n")
  cat("POST-REGISTRATION SNAPSHOT VALIDATION\n")
  cat("========================================\n")
  
  
  cat(
    "\nSnapshot folder exists: ",
    dir.exists(snapshot_root),
    "\n"
  )
  
  
  cat(
    "Final audit exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "00_admin/13_postregistration_exposure_replication_final_audit_v1.0.txt"
      )
    ),
    "\n"
  )
  
  
  cat(
    "Final R script exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "04_R_scripts/10_exploratory_exposure_replication.R"
      )
    ),
    "\n"
  )
  
  
  cat(
    "Final result directory exists: ",
    dir.exists(
      file.path(
        snapshot_root,
        "05_results/07_exposure_source_replication"
      )
    ),
    "\n"
  )
  
  
  cat(
    "Final figure PDF exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "06_figures/07_exposure_source_replication/Exposure_source_robustness_forest.pdf"
      )
    ),
    "\n"
  )
  
  
  cat(
    "Final figure PNG exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "06_figures/07_exposure_source_replication/Exposure_source_robustness_forest.png"
      )
    ),
    "\n"
  )
  
  
  cat(
    "File manifest exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "snapshot_file_manifest.csv"
      )
    ),
    "\n"
  )
  
  
  cat(
    "MD5 manifest exists: ",
    file.exists(
      file.path(
        snapshot_root,
        "snapshot_MD5_manifest.csv"
      )
    ),
    "\n"
  )
  
  
  cat(
    "\nTotal files in snapshot: ",
    length(
      list.files(
        snapshot_root,
        recursive = TRUE,
        full.names = TRUE
      )
    ),
    "\n"
  )
  
  
  cat(
    "\nSnapshot location:\n",
    snapshot_root,
    "\n"
  )