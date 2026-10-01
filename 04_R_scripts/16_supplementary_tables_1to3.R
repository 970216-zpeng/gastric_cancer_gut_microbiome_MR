# ============================================================
# Supplementary Table S1
# Eligibility assessment of candidate microbiome GWAS sources
# for the post-registration exposure-source robustness analyses
# ============================================================

library(here)
library(readr)
library(dplyr)
library(flextable)
library(officer)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- here(
  "07_tables",
  "02_publication_data",
  "TableS1_publication.csv"
)

output_file <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S1_external_GWAS_eligibility.docx"
)

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read publication data
# ------------------------------------------------------------

tab_s1 <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 3. Structural QC
# ------------------------------------------------------------

expected_cols <- c(
  "Data source",
  "Population",
  "Sample size",
  "Microbiome method",
  "Exact taxonomic match",
  "Full summary statistics",
  "Overlap concern",
  "Decision / rationale"
)

stopifnot(
  identical(names(tab_s1), expected_cols)
)

stopifnot(
  nrow(tab_s1) == 5
)

# Key inclusion/exclusion checks
stopifnot(
  grepl(
    "primary Post-reg",
    tab_s1$`Decision / rationale`[
      grepl("Swedish", tab_s1$`Data source`)
    ]
  )
)

stopifnot(
  grepl(
    "secondary Post-reg",
    tab_s1$`Decision / rationale`[
      grepl("Dutch Microbiome Project", tab_s1$`Data source`)
    ]
  )
)

stopifnot(
  grepl(
    "sample overlap",
    tab_s1$`Decision / rationale`[
      grepl("FINRISK", tab_s1$`Data source`)
    ]
  )
)

stopifnot(
  grepl(
    "no exact genus/family",
    tab_s1$`Decision / rationale`[
      grepl("HUNT", tab_s1$`Data source`)
    ]
  )
)

stopifnot(
  grepl(
    "cohort overlap",
    tab_s1$`Decision / rationale`[
      grepl("German", tab_s1$`Data source`)
    ]
  )
)

# ------------------------------------------------------------
# 4. Build flextable
# ------------------------------------------------------------

ft <- flextable(tab_s1)

ft <- font(
  ft,
  fontname = "Arial",
  part = "all"
)

ft <- fontsize(
  ft,
  size = 8,
  part = "all"
)

# Header
ft <- bold(
  ft,
  part = "header"
)

ft <- align(
  ft,
  align = "center",
  part = "header"
)

# Body
ft <- align(
  ft,
  j = c(
    "Data source",
    "Population",
    "Microbiome method",
    "Exact taxonomic match",
    "Full summary statistics",
    "Overlap concern",
    "Decision / rationale"
  ),
  align = "left",
  part = "body"
)

ft <- align(
  ft,
  j = "Sample size",
  align = "center",
  part = "body"
)

ft <- valign(
  ft,
  valign = "center",
  part = "all"
)

# ------------------------------------------------------------
# 5. Italicize Veillonella in Exact taxonomic match
# ------------------------------------------------------------

# Swedish
ft <- compose(
  ft,
  i = 1,
  j = "Exact taxonomic match",
  value = as_paragraph(
    as_chunk("Yes — exact genus "),
    as_i("Veillonella"),
    as_chunk(" and family Veillonellaceae")
  )
)

# DMP
ft <- compose(
  ft,
  i = 2,
  j = "Exact taxonomic match",
  value = as_paragraph(
    as_chunk("Yes — exact genus "),
    as_i("Veillonella"),
    as_chunk(" and family Veillonellaceae")
  )
)

# FINRISK
ft <- compose(
  ft,
  i = 3,
  j = "Exact taxonomic match",
  value = as_paragraph(
    as_chunk("Yes — exact "),
    as_i("Veillonella"),
    as_chunk(" and Veillonellaceae traits available in the frozen source audit")
  )
)

# HUNT
ft <- compose(
  ft,
  i = 4,
  j = "Exact taxonomic match",
  value = as_paragraph(
    as_chunk("No — available relevant "),
    as_i("Veillonella"),
    as_chunk(" traits were species-level; no exact genus "),
    as_i("Veillonella"),
    as_chunk(" / family Veillonellaceae pair identified in the frozen source audit")
  )
)

# ------------------------------------------------------------
# 6. Booktabs-style borders
# ------------------------------------------------------------

ft <- border_remove(ft)

border_top <- fp_border(
  color = "black",
  width = 1.2
)

border_mid <- fp_border(
  color = "black",
  width = 0.7
)

border_bottom <- fp_border(
  color = "black",
  width = 1.2
)

ft <- hline_top(
  ft,
  border = border_top,
  part = "header"
)

ft <- hline_bottom(
  ft,
  border = border_mid,
  part = "header"
)

ft <- hline_bottom(
  ft,
  border = border_bottom,
  part = "body"
)

# ------------------------------------------------------------
# 7. Padding
# ------------------------------------------------------------

ft <- padding(
  ft,
  padding.top = 3,
  padding.bottom = 3,
  padding.left = 3,
  padding.right = 3,
  part = "all"
)

# ------------------------------------------------------------
# 8. Column widths
# ------------------------------------------------------------

ft <- width(ft, j = "Data source",             width = 1.35)
ft <- width(ft, j = "Population",              width = 1.55)
ft <- width(ft, j = "Sample size",             width = 0.75)
ft <- width(ft, j = "Microbiome method",       width = 1.30)
ft <- width(ft, j = "Exact taxonomic match",   width = 2.15)
ft <- width(ft, j = "Full summary statistics", width = 2.10)
ft <- width(ft, j = "Overlap concern",         width = 2.25)
ft <- width(ft, j = "Decision / rationale",    width = 2.10)

# Keep rows intact
ft <- paginate(
  ft,
  init = TRUE,
  hdr_ftr = TRUE
)

# ------------------------------------------------------------
# 9. Title and abbreviation note
# ------------------------------------------------------------

title_text <- paste0(
  "Supplementary Table S1. Eligibility assessment of candidate ",
  "microbiome GWAS sources for the post-registration exposure-source ",
  "robustness analyses"
)

note1 <- paste0(
  "Abbreviations: GC, gastric cancer; GWAS, genome-wide association ",
  "study; Post-reg, post-registration."
)

note2 <- paste0(
  "Eligibility and exclusion decisions were finalized before inspection ",
  "of the corresponding external exposure-to-gastric-cancer MR results."
)

# ------------------------------------------------------------
# 10. Word document
# ------------------------------------------------------------

doc <- read_docx()

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      title_text,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc <- body_add_flextable(
  doc,
  value = ft
)

for (txt in c(note1, note2)) {
  
  doc <- body_add_fpar(
    doc,
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

doc <- body_end_section_landscape(doc)

print(
  doc,
  target = output_file
)

# ------------------------------------------------------------
# 11. Final QC report
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S1 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")
cat("Rows: ", nrow(tab_s1), "\n")
cat("Columns: ", ncol(tab_s1), "\n")
cat("Structural QC: PASS\n")
cat("Eligibility decision QC: PASS\n")
cat("Output:\n")
cat(output_file, "\n")
cat("====================================================\n")

# ============================================================
# SUPPLEMENTARY TABLE S2
# Phenotype definitions, processing, scales, and covariates
# ============================================================

# ------------------------------------------------------------
# S2.1 Paths
# ------------------------------------------------------------

s2_input <- here(
  "07_tables",
  "02_publication_data",
  "TableS2_publication.csv"
)

s2_output <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S2_phenotype_processing.docx"
)

# ------------------------------------------------------------
# S2.2 Read data
# ------------------------------------------------------------

tab_s2 <- read_csv(
  s2_input,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# S2.3 Structural QC
# ------------------------------------------------------------

expected_s2_cols <- c(
  "Data source",
  "GWAS ID / trait",
  "Phenotype definition",
  "Phenotype processing",
  "Association model / scale",
  "Covariate adjustment",
  "MR interpretation"
)

stopifnot(
  identical(names(tab_s2), expected_s2_cols)
)

stopifnot(
  nrow(tab_s2) == 6
)

# Dataset checks
stopifnot(
  any(grepl("MiBioGen", tab_s2$`Data source`))
)

stopifnot(
  any(grepl("FinnGen", tab_s2$`Data source`))
)

stopifnot(
  any(grepl("Swedish", tab_s2$`Data source`))
)

stopifnot(
  any(grepl("Dutch Microbiome Project", tab_s2$`Data source`))
)

# Key processing checks
stopifnot(
  grepl(
    "log",
    tab_s2$`Phenotype processing`[
      grepl("^MiBioGen$", tab_s2$`Data source`)
    ],
    ignore.case = TRUE
  )
)

stopifnot(
  grepl(
    "inverse-normal",
    tab_s2$`Phenotype processing`[
      grepl("Swedish", tab_s2$`Data source`)
    ],
    ignore.case = TRUE
  )
)

stopifnot(
  grepl(
    "inverse-rank",
    tab_s2$`Phenotype processing`[
      grepl("Dutch Microbiome Project", tab_s2$`Data source`)
    ],
    ignore.case = TRUE
  )
)

# ------------------------------------------------------------
# S2.4 Build flextable
# ------------------------------------------------------------

ft_s2 <- flextable(tab_s2)

ft_s2 <- font(
  ft_s2,
  fontname = "Arial",
  part = "all"
)

ft_s2 <- fontsize(
  ft_s2,
  size = 7.5,
  part = "all"
)

ft_s2 <- bold(
  ft_s2,
  part = "header"
)

ft_s2 <- align(
  ft_s2,
  align = "center",
  part = "header"
)

ft_s2 <- align(
  ft_s2,
  align = "left",
  part = "body"
)

ft_s2 <- valign(
  ft_s2,
  valign = "center",
  part = "all"
)

# ------------------------------------------------------------
# S2.5 Booktabs borders
# ------------------------------------------------------------

ft_s2 <- border_remove(ft_s2)

s2_border_top <- fp_border(
  color = "black",
  width = 1.2
)

s2_border_mid <- fp_border(
  color = "black",
  width = 0.7
)

s2_border_bottom <- fp_border(
  color = "black",
  width = 1.2
)

ft_s2 <- hline_top(
  ft_s2,
  border = s2_border_top,
  part = "header"
)

ft_s2 <- hline_bottom(
  ft_s2,
  border = s2_border_mid,
  part = "header"
)

ft_s2 <- hline_bottom(
  ft_s2,
  border = s2_border_bottom,
  part = "body"
)

# ------------------------------------------------------------
# S2.6 Padding
# ------------------------------------------------------------

ft_s2 <- padding(
  ft_s2,
  padding.top = 3,
  padding.bottom = 3,
  padding.left = 3,
  padding.right = 3,
  part = "all"
)

# ------------------------------------------------------------
# S2.7 Column widths
# ------------------------------------------------------------

ft_s2 <- width(ft_s2, j = "Data source",                width = 1.45)
ft_s2 <- width(ft_s2, j = "GWAS ID / trait",            width = 1.70)
ft_s2 <- width(ft_s2, j = "Phenotype definition",       width = 2.10)
ft_s2 <- width(ft_s2, j = "Phenotype processing",       width = 2.55)
ft_s2 <- width(ft_s2, j = "Association model / scale",  width = 2.20)
ft_s2 <- width(ft_s2, j = "Covariate adjustment",       width = 2.55)
ft_s2 <- width(ft_s2, j = "MR interpretation",          width = 2.35)

ft_s2 <- paginate(
  ft_s2,
  init = TRUE,
  hdr_ftr = TRUE
)

# ------------------------------------------------------------
# S2.8 Title and notes
# ------------------------------------------------------------

s2_title <- paste0(
  "Supplementary Table S2. Phenotype definitions, source-specific ",
  "processing, association scales, and covariate adjustment of the ",
  "GWAS datasets used in the Mendelian randomization analyses"
)

s2_note1 <- paste0(
  "Abbreviations: GC, gastric cancer; GWAS, genome-wide association study; ",
  "MR, Mendelian randomization; PC, principal component; ",
  "RIN, rank-based inverse-normal transformation."
)

s2_note2 <- paste0(
  "Microbiome GWAS effect estimates reflect source-specific transformed ",
  "abundance phenotypes. Effect units across MiBioGen, Swedish 2026, and ",
  "DMP should therefore not be assumed to be directly commensurate."
)

# ------------------------------------------------------------
# S2.9 Word output
# ------------------------------------------------------------

doc_s2 <- read_docx()

doc_s2 <- body_add_fpar(
  doc_s2,
  fpar(
    ftext(
      s2_title,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc_s2 <- body_add_flextable(
  doc_s2,
  value = ft_s2
)

for (txt in c(s2_note1, s2_note2)) {
  
  doc_s2 <- body_add_fpar(
    doc_s2,
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

doc_s2 <- body_end_section_landscape(doc_s2)

print(
  doc_s2,
  target = s2_output
)

# ------------------------------------------------------------
# S2.10 QC report
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S2 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")
cat("Rows: ", nrow(tab_s2), "\n")
cat("Columns: ", ncol(tab_s2), "\n")
cat("Structural QC: PASS\n")
cat("Phenotype-processing QC: PASS\n")
cat("Output:\n")
cat(s2_output, "\n")
cat("====================================================\n")

# ============================================================
# SUPPLEMENTARY TABLE S3
# Source-file discovery for instrument flow and SNP-level data
# ============================================================

library(here)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)

# ------------------------------------------------------------
# S3.A1 Search roots
# ------------------------------------------------------------

s3_search_roots <- c(
  here("02_intermediate_data"),
  here("03_clean_data"),
  here("05_results"),
  here("09_snapshot")
)

s3_search_roots <- s3_search_roots[
  dir.exists(s3_search_roots)
]

stopifnot(length(s3_search_roots) > 0)

# ------------------------------------------------------------
# S3.A2 Find candidate files
# ------------------------------------------------------------

s3_files <- unlist(
  lapply(
    s3_search_roots,
    function(x) {
      list.files(
        x,
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE
      )
    }
  )
)

s3_files <- unique(s3_files)

# Keep likely analytical files only
s3_files <- s3_files[
  grepl(
    "\\.(csv|tsv|txt|rds)$",
    s3_files,
    ignore.case = TRUE
  )
]

# Filename relevance
s3_name_flag <- grepl(
  paste0(
    "harm|single|snp|instrument|clump|",
    "forward|reverse|veillon|gcst|",
    "primary|threshold|finn|east|swedish|dmp"
  ),
  basename(s3_files),
  ignore.case = TRUE
)

s3_files <- s3_files[s3_name_flag]

# ------------------------------------------------------------
# S3.A3 Safely inspect column names
# ------------------------------------------------------------

get_columns_safe <- function(path) {
  
  ext <- tolower(tools::file_ext(path))
  
  out <- tryCatch({
    
    if (ext == "csv") {
      
      x <- suppressMessages(
        read_csv(
          path,
          n_max = 0,
          show_col_types = FALSE,
          progress = FALSE
        )
      )
      
      names(x)
      
    } else if (ext %in% c("tsv", "txt")) {
      
      x <- suppressMessages(
        read_tsv(
          path,
          n_max = 0,
          show_col_types = FALSE,
          progress = FALSE
        )
      )
      
      names(x)
      
    } else if (ext == "rds") {
      
      obj <- readRDS(path)
      
      if (is.data.frame(obj)) {
        names(obj)
      } else {
        character(0)
      }
      
    } else {
      
      character(0)
    }
    
  }, error = function(e) {
    
    character(0)
    
  })
  
  out
}

# ------------------------------------------------------------
# S3.A4 Build inventory
# ------------------------------------------------------------

s3_inventory <- map_dfr(
  s3_files,
  function(path) {
    
    cols <- get_columns_safe(path)
    
    cols_lower <- tolower(cols)
    
    tibble(
      file_name = basename(path),
      relative_path = str_replace(
        normalizePath(path, winslash = "/", mustWork = FALSE),
        paste0(
          "^",
          str_replace_all(
            normalizePath(here(), winslash = "/", mustWork = FALSE),
            "([\\W])",
            "\\\\\\1"
          ),
          "/?"
        ),
        ""
      ),
      
      n_columns = length(cols),
      
      has_SNP =
        any(cols_lower %in% c("snp", "rsid", "variant", "variant_id")),
      
      has_beta_exposure =
        any(grepl(
          "beta.*exposure|exposure.*beta",
          cols_lower
        )),
      
      has_se_exposure =
        any(grepl(
          "se.*exposure|exposure.*se",
          cols_lower
        )),
      
      has_beta_outcome =
        any(grepl(
          "beta.*outcome|outcome.*beta",
          cols_lower
        )),
      
      has_se_outcome =
        any(grepl(
          "se.*outcome|outcome.*se",
          cols_lower
        )),
      
      has_pval_exposure =
        any(grepl(
          "pval.*exposure|p.*exposure",
          cols_lower
        )),
      
      has_mr_keep =
        any(cols_lower == "mr_keep"),
      
      has_effect_allele =
        any(grepl(
          "effect_allele",
          cols_lower
        )),
      
      has_other_allele =
        any(grepl(
          "other_allele",
          cols_lower
        )),
      
      columns = paste(cols, collapse = " | ")
    )
  }
)

# ------------------------------------------------------------
# S3.A5 Candidate scoring
# ------------------------------------------------------------

s3_inventory <- s3_inventory %>%
  mutate(
    snp_level_score =
      has_SNP +
      has_beta_exposure +
      has_se_exposure +
      has_beta_outcome +
      has_se_outcome +
      has_mr_keep,
    
    likely_harmonised =
      has_SNP &
      has_beta_exposure &
      has_se_exposure &
      has_beta_outcome &
      has_se_outcome,
    
    likely_final_MR =
      likely_harmonised &
      has_mr_keep
  ) %>%
  arrange(
    desc(likely_final_MR),
    desc(likely_harmonised),
    desc(snp_level_score),
    file_name
  )

# ------------------------------------------------------------
# S3.A6 Save audit inventory
# ------------------------------------------------------------

s3_inventory_output <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_source_file_inventory.csv"
)

dir.create(
  dirname(s3_inventory_output),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  s3_inventory,
  s3_inventory_output
)

# ------------------------------------------------------------
# S3.A7 Console report
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TABLE S3 SOURCE-FILE DISCOVERY COMPLETE\n")
cat("====================================================\n")
cat("Candidate files: ", nrow(s3_inventory), "\n")
cat(
  "Likely harmonised SNP files: ",
  sum(s3_inventory$likely_harmonised),
  "\n"
)
cat(
  "Likely final MR files: ",
  sum(s3_inventory$likely_final_MR),
  "\n"
)
cat("Inventory:\n")
cat(s3_inventory_output, "\n")
cat("====================================================\n")

print(
  s3_inventory %>%
    filter(
      likely_harmonised |
        snp_level_score >= 3
    ) %>%
    select(
      file_name,
      relative_path,
      snp_level_score,
      likely_harmonised,
      likely_final_MR
    ),
  n = Inf
)
# ============================================================
# SUPPLEMENTARY TABLE S3
# Canonical source loading, re-harmonisation validation,
# instrument flow, and SNP-level audit master
# ============================================================

library(here)
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(TwoSampleMR)

# ------------------------------------------------------------
# S3.B1 Source map
# ------------------------------------------------------------

s3_map_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_analysis_source_map.csv"
)

s3_map <- read_csv(
  s3_map_file,
  show_col_types = FALSE
)

stopifnot(
  nrow(s3_map) == 18,
  !anyDuplicated(s3_map$analysis_id)
)

# ------------------------------------------------------------
# S3.B2 Helper: read CSV / RDS
# ------------------------------------------------------------

read_s3_file <- function(rel_path) {
  
  full_path <- here(rel_path)
  
  if (!file.exists(full_path)) {
    stop("Missing source file: ", rel_path)
  }
  
  ext <- tolower(tools::file_ext(full_path))
  
  if (ext == "csv") {
    
    read_csv(
      full_path,
      show_col_types = FALSE
    )
    
  } else if (ext == "rds") {
    
    readRDS(full_path)
    
  } else {
    
    stop("Unsupported file type: ", full_path)
  }
}

# ------------------------------------------------------------
# S3.B3 Convert generic extracted outcome data to
#       TwoSampleMR outcome format
# ------------------------------------------------------------

format_s3_outcome <- function(
    raw,
    outcome_name,
    outcome_id
) {
  
  required <- c(
    "SNP",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "pval"
  )
  
  missing_cols <- setdiff(
    required,
    names(raw)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing outcome columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- tibble(
    SNP = as.character(raw$SNP),
    
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
    
    outcome = outcome_name,
    
    id.outcome = outcome_id,
    
    mr_keep.outcome = TRUE,
    
    pval_origin.outcome = "reported"
  )
  
  if ("eaf" %in% names(raw)) {
    out$eaf.outcome <- as.numeric(raw$eaf)
  }
  
  if ("samplesize" %in% names(raw)) {
    out$samplesize.outcome <- as.numeric(raw$samplesize)
  }
  
  if ("ncase" %in% names(raw)) {
    out$ncase.outcome <- as.numeric(raw$ncase)
  }
  
  if ("chr" %in% names(raw)) {
    out$chr.outcome <- raw$chr
  }
  
  if ("pos" %in% names(raw)) {
    out$pos.outcome <- raw$pos
  }
  
  out
}

# ------------------------------------------------------------
# S3.B4 Locked expected final SNP counts
# ------------------------------------------------------------

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
    s3_map$analysis_id
  )
)
# ------------------------------------------------------------
# Locked external pre-harmonisation instrument counts
# ------------------------------------------------------------

expected_external_selected_n <- c(
  
  swedish_primary_V = 7,
  swedish_primary_F = 8,
  
  dmp_primary_V = 5,
  dmp_primary_F = 7,
  
  swedish_p1e5_V = 15,
  swedish_p1e5_F = 18,
  
  dmp_p1e5_V = 8,
  dmp_p1e5_F = 13
)

# ------------------------------------------------------------
# S3.B5 Process one analysis
# ------------------------------------------------------------

process_s3_analysis <- function(i) {
  
  meta <- s3_map[i, ]
  
  mode <- meta$source_mode
  
  # ----------------------------------------------------------
  # A. Direct harmonised data
  # ----------------------------------------------------------
  
  if (mode == "direct_harmonised") {
    
    dat <- read_s3_file(
      meta$harmonised_path
    )
    
    # Number entering the harmonisation/outcome stage
    available_n <- nrow(dat)
    
    # Where the frozen exposure IV file is available,
    # use it to count the post-clumping/F>10 instrument set.
    if (
      !is.na(meta$exposure_path) &&
      nzchar(meta$exposure_path)
    ) {
      
      # Core analyses: count frozen clumped exposure IVs directly
      exposure_dat <- read_s3_file(
        meta$exposure_path
      )
      
      selected_n <- nrow(exposure_dat)
      
      # For these core direct-harmonised datasets,
      # available SNP count is represented by the harmonised input.
      available_n <- nrow(dat)
      
    } else if (
      meta$analysis_id %in%
      names(expected_external_selected_n)
    ) {
      
      # Post-registration external analyses:
      # frozen final audit supplies the pre-harmonisation clumped/F>10 count.
      selected_n <- unname(
        expected_external_selected_n[
          meta$analysis_id
        ]
      )
      
      # Frozen outcome-extraction audit:
      # all requested external IV SNPs were found.
      available_n <- selected_n
      
    } else {
      
      stop(
        "No canonical pre-harmonisation count source for: ",
        meta$analysis_id
      )
    }
    
    # ----------------------------------------------------------
    # B. Re-harmonise frozen exposure + extracted outcome
    # ----------------------------------------------------------
    
  } else if (
    mode %in% c(
      "reharmonise",
      "reharmonise_reverse"
    )
  ) {
    
    exposure_dat <- read_s3_file(
      meta$exposure_path
    )
    
    outcome_raw <- read_s3_file(
      meta$outcome_path
    )
    
    selected_n <- nrow(exposure_dat)
    
    # Label outcome according to analysis direction
    if (mode == "reharmonise_reverse") {
      
      out_name <- meta$outcome
      out_id <- meta$outcome
      
    } else {
      
      out_name <- meta$outcome
      out_id <- meta$outcome
    }
    
    outcome_dat <- format_s3_outcome(
      raw = outcome_raw,
      outcome_name = out_name,
      outcome_id = out_id
    )
    
    # Outcome availability before harmonisation
    available_n <- sum(
      exposure_dat$SNP %in%
        outcome_dat$SNP
    )
    
    # Reproduce frozen harmonisation rule
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = 3
    )
    
  } else {
    
    stop(
      "Unknown source_mode: ",
      mode
    )
  }
  
  # ----------------------------------------------------------
  # C. Structural checks
  # ----------------------------------------------------------
  
  required_harmonised <- c(
    "SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "beta.exposure",
    "se.exposure",
    "beta.outcome",
    "se.outcome",
    "pval.exposure",
    "pval.outcome",
    "mr_keep"
  )
  
  missing_harm_cols <- setdiff(
    required_harmonised,
    names(dat)
  )
  
  if (length(missing_harm_cols) > 0) {
    
    stop(
      meta$analysis_id,
      ": missing harmonised columns: ",
      paste(
        missing_harm_cols,
        collapse = ", "
      )
    )
  }
  
  # normalize logical field
  dat$mr_keep <- as.logical(
    dat$mr_keep
  )
  
  # Instrument strength
  dat <- dat %>%
    mutate(
      F_stat_S3 =
        (beta.exposure /
           se.exposure)^2
    )
  
  final_dat <- dat %>%
    filter(mr_keep %in% TRUE)
  
  final_n <- nrow(final_dat)
  
  # ----------------------------------------------------------
  # D. Locked-count validation
  # ----------------------------------------------------------
  
  expected_n <-
    unname(
      expected_final_n[
        meta$analysis_id
      ]
    )
  
  if (final_n != expected_n) {
    
    stop(
      "\nS3 FINAL SNP COUNT MISMATCH\n",
      "Analysis: ",
      meta$analysis_id,
      "\nExpected: ",
      expected_n,
      "\nObserved: ",
      final_n,
      "\n"
    )
  }
  
  if (anyDuplicated(final_dat$SNP)) {
    
    stop(
      "Duplicate final SNPs in: ",
      meta$analysis_id
    )
  }
  
  # ----------------------------------------------------------
  # E. Instrument-flow row
  # ----------------------------------------------------------
  
  flow <- tibble(
    analysis_id =
      meta$analysis_id,
    
    analysis_group =
      meta$analysis_group,
    
    direction =
      meta$direction,
    
    exposure_source =
      meta$exposure_source,
    
    exposure =
      meta$exposure,
    
    outcome =
      meta$outcome,
    
    instrument_threshold =
      meta$instrument_threshold,
    
    source_mode =
      meta$source_mode,
    
    selected_clumped_Fgt10 =
      selected_n,
    
    available_in_outcome =
      available_n,
    
    missing_in_outcome =
      selected_n -
      available_n,
    
    harmonised_rows =
      nrow(dat),
    
    final_MR_SNPs =
      final_n,
    
    excluded_after_availability =
      available_n -
      final_n,
    
    expected_final_MR_SNPs =
      expected_n,
    
    final_count_check =
      "PASS"
  )
  
  # ----------------------------------------------------------
  # F. Final SNP-level rows
  # ----------------------------------------------------------
  
  snp <- final_dat %>%
    transmute(
      
      analysis_id =
        meta$analysis_id,
      
      analysis_group =
        meta$analysis_group,
      
      direction =
        meta$direction,
      
      exposure_source =
        meta$exposure_source,
      
      exposure =
        meta$exposure,
      
      outcome =
        meta$outcome,
      
      instrument_threshold =
        meta$instrument_threshold,
      
      SNP =
        SNP,
      
      effect_allele =
        effect_allele.exposure,
      
      other_allele =
        other_allele.exposure,
      
      beta_exposure =
        beta.exposure,
      
      se_exposure =
        se.exposure,
      
      p_exposure =
        pval.exposure,
      
      beta_outcome =
        beta.outcome,
      
      se_outcome =
        se.outcome,
      
      p_outcome =
        pval.outcome,
      
      F_stat =
        F_stat_S3,
      
      mr_keep =
        mr_keep
    )
  
  list(
    flow = flow,
    snp = snp,
    harmonised = final_dat
  )
}

# ------------------------------------------------------------
# S3.B6 Run all 18 analyses
# ------------------------------------------------------------

s3_processed <- map(
  seq_len(nrow(s3_map)),
  process_s3_analysis
)

s3_flow <- bind_rows(
  map(
    s3_processed,
    "flow"
  )
)

s3_snp <- bind_rows(
  map(
    s3_processed,
    "snp"
  )
)

# ------------------------------------------------------------
# S3.B7 Global QC
# ------------------------------------------------------------

stopifnot(
  nrow(s3_flow) == 18
)

stopifnot(
  all(
    s3_flow$final_count_check ==
      "PASS"
  )
)

stopifnot(
  all(
    s3_snp$F_stat > 10
  )
)

stopifnot(
  all(
    s3_snp$mr_keep
  )
)

# SNP count must equal flow totals
s3_snp_counts <- s3_snp %>%
  count(
    analysis_id,
    name = "n_snp_rows"
  )

s3_flow_check <- s3_flow %>%
  left_join(
    s3_snp_counts,
    by = "analysis_id"
  )

stopifnot(
  all(
    s3_flow_check$final_MR_SNPs ==
      s3_flow_check$n_snp_rows
  )
)

# ------------------------------------------------------------
# S3.B8 Re-harmonisation-specific validation
# ------------------------------------------------------------

s3_reharmonised_qc <- s3_flow %>%
  filter(
    source_mode !=
      "direct_harmonised"
  ) %>%
  select(
    analysis_id,
    selected_clumped_Fgt10,
    available_in_outcome,
    harmonised_rows,
    final_MR_SNPs,
    expected_final_MR_SNPs,
    final_count_check
  )

stopifnot(
  nrow(s3_reharmonised_qc) == 6
)

# ------------------------------------------------------------
# S3.B9 Save audit outputs
# ------------------------------------------------------------

s3_flow_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_instrument_flow_complete_audit.csv"
)

s3_snp_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_SNP_level_audit_master.csv"
)

s3_reharm_file <- here(
  "07_tables",
  "01_audit_data",
  "TableS3_reharmonisation_validation.csv"
)

write_csv(
  s3_flow,
  s3_flow_file
)

write_csv(
  s3_snp,
  s3_snp_file
)

write_csv(
  s3_reharmonised_qc,
  s3_reharm_file
)

# ------------------------------------------------------------
# S3.B10 Console QC report
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S3 DATA AUDIT COMPLETE\n")
cat("====================================================\n")

cat(
  "Analyses processed: ",
  nrow(s3_flow),
  "\n"
)

cat(
  "Direct harmonised analyses: ",
  sum(
    s3_flow$source_mode ==
      "direct_harmonised"
  ),
  "\n"
)

cat(
  "Re-harmonised analyses: ",
  sum(
    s3_flow$source_mode !=
      "direct_harmonised"
  ),
  "\n"
)

cat(
  "Total final SNP-analysis rows: ",
  nrow(s3_snp),
  "\n"
)

cat(
  "Final SNP count checks: PASS\n"
)

cat(
  "F-statistic > 10 checks: PASS\n"
)

cat(
  "Duplicate-within-analysis checks: PASS\n"
)

cat("\nInstrument-flow audit:\n")
cat(
  s3_flow_file,
  "\n"
)

cat("\nSNP-level audit master:\n")
cat(
  s3_snp_file,
  "\n"
)

cat("\nRe-harmonisation validation:\n")
cat(
  s3_reharm_file,
  "\n"
)

cat("====================================================\n")

print(
  s3_flow %>%
    select(
      analysis_id,
      selected_clumped_Fgt10,
      available_in_outcome,
      harmonised_rows,
      final_MR_SNPs,
      final_count_check
    ),
  n = Inf
)
# ============================================================
# SUPPLEMENTARY TABLE S3
# Publication layer + Word output
# ============================================================

library(dplyr)
library(readr)
library(flextable)
library(officer)
library(stringr)

# ------------------------------------------------------------
# S3.C1 Analysis labels
# ------------------------------------------------------------

s3_labels <- c(
  
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
    names(s3_labels),
    s3_flow$analysis_id
  )
)

# ============================================================
# PANEL A — Instrument flow
# ============================================================

s3A_pub <- s3_flow %>%
  mutate(
    Analysis = unname(
      s3_labels[analysis_id]
    ),
    
    Threshold = instrument_threshold
  ) %>%
  transmute(
    Analysis,
    Threshold,
    `Clumped / F > 10` =
      selected_clumped_Fgt10,
    
    `Available in outcome` =
      available_in_outcome,
    
    `Harmonised rows` =
      harmonised_rows,
    
    `Final MR SNPs` =
      final_MR_SNPs
  )

# QC
stopifnot(
  nrow(s3A_pub) == 18
)

# Save publication CSV
s3A_csv <- here(
  "07_tables",
  "02_publication_data",
  "TableS3A_instrument_flow_publication.csv"
)

write_csv(
  s3A_pub,
  s3A_csv
)

# ============================================================
# PANEL B — SNP-level associations
# ============================================================

format_beta_se <- function(beta, se) {
  
  sprintf(
    "%.5f (%.5f)",
    beta,
    se
  )
}

format_p <- function(p) {
  
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
      sprintf("%.4f", p)
    )
  )
}

s3B_pub <- s3_snp %>%
  mutate(
    
    Analysis = unname(
      s3_labels[analysis_id]
    ),
    
    Threshold =
      instrument_threshold,
    
    `EA / OA` =
      paste0(
        effect_allele,
        " / ",
        other_allele
      ),
    
    `β exposure (SE)` =
      format_beta_se(
        beta_exposure,
        se_exposure
      ),
    
    `P exposure` =
      format_p(
        p_exposure
      ),
    
    `β outcome (SE)` =
      format_beta_se(
        beta_outcome,
        se_outcome
      ),
    
    `P outcome` =
      format_p(
        p_outcome
      ),
    
    `F statistic` =
      sprintf(
        "%.1f",
        F_stat
      )
  ) %>%
  transmute(
    Analysis,
    Threshold,
    SNP,
    `EA / OA`,
    `β exposure (SE)`,
    `P exposure`,
    `β outcome (SE)`,
    `P outcome`,
    `F statistic`
  )

# ------------------------------------------------------------
# S3.C2 SNP-level QC
# ------------------------------------------------------------

stopifnot(
  nrow(s3B_pub) ==
    sum(s3_flow$final_MR_SNPs)
)

stopifnot(
  nrow(s3B_pub) == 125
)

stopifnot(
  !anyNA(s3B_pub$SNP)
)

stopifnot(
  !anyNA(s3B_pub$`F statistic`)
)

# Save publication CSV
s3B_csv <- here(
  "07_tables",
  "02_publication_data",
  "TableS3B_SNP_level_publication.csv"
)

write_csv(
  s3B_pub,
  s3B_csv
)

# ============================================================
# S3.C3 PANEL A flextable
# ============================================================

ft_s3A <- flextable(
  s3A_pub
)

ft_s3A <- font(
  ft_s3A,
  fontname = "Arial",
  part = "all"
)

ft_s3A <- fontsize(
  ft_s3A,
  size = 8,
  part = "all"
)

ft_s3A <- bold(
  ft_s3A,
  part = "header"
)

ft_s3A <- align(
  ft_s3A,
  align = "center",
  part = "header"
)

ft_s3A <- align(
  ft_s3A,
  j = "Analysis",
  align = "left",
  part = "body"
)

ft_s3A <- align(
  ft_s3A,
  j = 2:ncol(s3A_pub),
  align = "center",
  part = "body"
)

ft_s3A <- valign(
  ft_s3A,
  valign = "center",
  part = "all"
)

# Booktabs
ft_s3A <- border_remove(
  ft_s3A
)

ft_s3A <- hline_top(
  ft_s3A,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "header"
)

ft_s3A <- hline_bottom(
  ft_s3A,
  border = fp_border(
    color = "black",
    width = 0.7
  ),
  part = "header"
)

ft_s3A <- hline_bottom(
  ft_s3A,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "body"
)

ft_s3A <- padding(
  ft_s3A,
  padding.top = 2.5,
  padding.bottom = 2.5,
  padding.left = 3,
  padding.right = 3,
  part = "all"
)

ft_s3A <- width(
  ft_s3A,
  j = "Analysis",
  width = 4.0
)

ft_s3A <- width(
  ft_s3A,
  j = "Threshold",
  width = 1.35
)

ft_s3A <- width(
  ft_s3A,
  j = "Clumped / F > 10",
  width = 1.25
)

ft_s3A <- width(
  ft_s3A,
  j = "Available in outcome",
  width = 1.40
)

ft_s3A <- width(
  ft_s3A,
  j = "Harmonised rows",
  width = 1.25
)

ft_s3A <- width(
  ft_s3A,
  j = "Final MR SNPs",
  width = 1.15
)

ft_s3A <- paginate(
  ft_s3A,
  init = TRUE,
  hdr_ftr = TRUE
)

# ============================================================
# S3.C4 PANEL B flextable
# ============================================================

ft_s3B <- flextable(
  s3B_pub
)

ft_s3B <- font(
  ft_s3B,
  fontname = "Arial",
  part = "all"
)

ft_s3B <- fontsize(
  ft_s3B,
  size = 7,
  part = "all"
)

ft_s3B <- bold(
  ft_s3B,
  part = "header"
)

ft_s3B <- align(
  ft_s3B,
  align = "center",
  part = "header"
)

ft_s3B <- align(
  ft_s3B,
  j = "Analysis",
  align = "left",
  part = "body"
)

ft_s3B <- align(
  ft_s3B,
  j = 2:ncol(s3B_pub),
  align = "center",
  part = "body"
)

ft_s3B <- valign(
  ft_s3B,
  valign = "center",
  part = "all"
)

# Borders
ft_s3B <- border_remove(
  ft_s3B
)

ft_s3B <- hline_top(
  ft_s3B,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "header"
)

ft_s3B <- hline_bottom(
  ft_s3B,
  border = fp_border(
    color = "black",
    width = 0.7
  ),
  part = "header"
)

ft_s3B <- hline_bottom(
  ft_s3B,
  border = fp_border(
    color = "black",
    width = 1.2
  ),
  part = "body"
)

ft_s3B <- padding(
  ft_s3B,
  padding.top = 1.7,
  padding.bottom = 1.7,
  padding.left = 2,
  padding.right = 2,
  part = "all"
)

# Widths
ft_s3B <- width(
  ft_s3B,
  j = "Analysis",
  width = 3.3
)

ft_s3B <- width(
  ft_s3B,
  j = "Threshold",
  width = 1.10
)

ft_s3B <- width(
  ft_s3B,
  j = "SNP",
  width = 1.05
)

ft_s3B <- width(
  ft_s3B,
  j = "EA / OA",
  width = 0.70
)

ft_s3B <- width(
  ft_s3B,
  j = "β exposure (SE)",
  width = 1.35
)

ft_s3B <- width(
  ft_s3B,
  j = "P exposure",
  width = 1.05
)

ft_s3B <- width(
  ft_s3B,
  j = "β outcome (SE)",
  width = 1.35
)

ft_s3B <- width(
  ft_s3B,
  j = "P outcome",
  width = 1.05
)

ft_s3B <- width(
  ft_s3B,
  j = "F statistic",
  width = 0.80
)

ft_s3B <- set_table_properties(
  ft_s3B,
  opts_word = list(
    split = TRUE,
    repeat_headers = TRUE
  )
)

# ============================================================
# S3.C5 Word output
# ============================================================

s3_output <- here(
  "07_tables",
  "03_outputs",
  "Supplementary_Table_S3_instrument_flow_SNP_level.docx"
)

s3_title <- paste0(
  "Supplementary Table S3. Instrument flow and SNP-level ",
  "genetic associations for the Mendelian randomization analyses"
)

s3_note1 <- paste0(
  "Abbreviations: Alt-E, alternative European; DMP, Dutch Microbiome ",
  "Project; EA, effect allele; GC, gastric cancer; MR, Mendelian ",
  "randomization; OA, other allele; Post-reg, post-registration; ",
  "SE, standard error; SNP, single-nucleotide polymorphism."
)

s3_note2 <- paste0(
  "All instruments satisfied F > 10. No proxy variants were used. ",
  "Exposure and outcome alleles were harmonized using the conservative ",
  "TwoSampleMR action = 3 procedure."
)

s3_note3 <- paste0(
  "Panel A shows instrument flow from the clumped, strength-filtered ",
  "instrument set through outcome availability and harmonization to ",
  "the final MR instrument set. Panel B reports the final harmonized ",
  "SNP-level associations entering MR."
)

doc_s3 <- read_docx()

doc_s3 <- body_add_fpar(
  doc_s3,
  fpar(
    ftext(
      s3_title,
      prop = fp_text(
        font.family = "Arial",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

# Panel A
doc_s3 <- body_add_fpar(
  doc_s3,
  fpar(
    ftext(
      "Panel A. Instrument flow",
      prop = fp_text(
        font.family = "Arial",
        font.size = 9,
        bold = TRUE
      )
    )
  )
)

doc_s3 <- body_add_flextable(
  doc_s3,
  value = ft_s3A
)

# Page break before long SNP table
doc_s3 <- body_add_break(
  doc_s3
)

# Panel B
doc_s3 <- body_add_fpar(
  doc_s3,
  fpar(
    ftext(
      "Panel B. SNP-level genetic associations",
      prop = fp_text(
        font.family = "Arial",
        font.size = 9,
        bold = TRUE
      )
    ),
    fp_p = fp_par(
      keep_with_next = TRUE
    )
  )
)

doc_s3 <- body_add_flextable(
  doc_s3,
  value = ft_s3B
)

# Notes
for (txt in c(
  s3_note1,
  s3_note2,
  s3_note3
)) {
  
  doc_s3 <- body_add_fpar(
    doc_s3,
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

doc_s3 <- body_end_section_landscape(
  doc_s3
)

print(
  doc_s3,
  target = s3_output
)

# ============================================================
# S3.C6 Final report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("SUPPLEMENTARY TABLE S3 GENERATED SUCCESSFULLY\n")
cat("====================================================\n")
cat("Panel A rows: ", nrow(s3A_pub), "\n")
cat("Panel B rows: ", nrow(s3B_pub), "\n")
cat("Expected SNP rows: 125\n")
cat("Instrument-flow QC: PASS\n")
cat("SNP-level QC: PASS\n")
cat("Final-count reconciliation: PASS\n")
cat("Output:\n")
cat(s3_output, "\n")
cat("====================================================\n")
