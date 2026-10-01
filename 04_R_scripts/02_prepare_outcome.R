# ============================================================
# 02_prepare_outcome.R
# Project: 01_MR_GastricCancer
# Primary outcome: Gastric cancer
# GWAS: GCST90018849
# ============================================================

# 1. Load packages --------------------------------------------

library(TwoSampleMR)


# 2. Load primary exposure IVs --------------------------------

veillonella_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

veillonellaceae_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_primary_p5e6.rds"
)


# Check
nrow(veillonella_iv)
nrow(veillonellaceae_iv)
# 3. Create SNP list ------------------------------------------

target_snps <- unique(
  c(
    veillonella_iv$SNP,
    veillonellaceae_iv$SNP
  )
)

length(target_snps)

target_snps
# 4. Define outcome VCF ---------------------------------------

file_gc_primary <-
  "01_raw_data/02_outcome/ebi-a-GCST90018849.vcf.gz"

file.exists(file_gc_primary)
# 5. Inspect gastric cancer VCF -------------------------------

con <- gzfile(file_gc_primary, open = "rt")

gc_head <- readLines(con, n = 200)

close(con)


# Column header
gc_head[grepl("^#CHROM", gc_head)]

# FORMAT definitions
gc_head[grepl("^##FORMAT", gc_head)]
# ============================================================
# 6. Extract target IVs from primary gastric cancer outcome
# ============================================================

gc_primary_result <- extract_outcome_snps(
  vcf_file = file_gc_primary,
  target_snps = target_snps,
  chunk_size = 50000
)

# Extract data object
gc_primary <- gc_primary_result$data
# ============================================================
# 7. Check SNP matching
# ============================================================

gc_primary_result$n_requested
gc_primary_result$n_found
gc_primary_result$n_missing

gc_primary_result$missing
# ============================================================
# 8. Inspect extracted outcome associations
# ============================================================

gc_primary[
  ,
  c(
    "SNP",
    "chr",
    "pos",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "pval",
    "eaf",
    "samplesize",
    "ncase"
  )
]
# ============================================================
# 9. Check data completeness
# ============================================================

sum(is.na(gc_primary$beta))
sum(is.na(gc_primary$se))
sum(is.na(gc_primary$pval))
sum(is.na(gc_primary$eaf))
sum(is.na(gc_primary$samplesize))
sum(is.na(gc_primary$ncase))

sum(duplicated(gc_primary$SNP))
# ============================================================
# 10. Save primary outcome data
# ============================================================

dir.create(
  "02_intermediate_data/02_outcome",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  gc_primary,
  "02_intermediate_data/02_outcome/GCST90018849_primary_outcome_IVs.csv",
  row.names = FALSE
)

saveRDS(
  gc_primary,
  "02_intermediate_data/02_outcome/GCST90018849_primary_outcome_IVs.rds"
)

#检查暴露
names(veillonella_iv)

names(veillonellaceae_iv)


"eaf.exposure" %in% names(veillonella_iv)

"eaf.exposure" %in% names(veillonellaceae_iv)
# ============================================================
# 11. Check exposure EAF completeness
# ============================================================

sum(is.na(veillonella_iv$eaf.exposure))
sum(is.na(veillonellaceae_iv$eaf.exposure))

veillonella_iv[, c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "eaf.exposure"
)]

veillonellaceae_iv[, c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "eaf.exposure"
)]
# ============================================================
# 12. Recover exposure EAF from original GWAS-VCF
# ============================================================

veillonella_eaf_result <- extract_outcome_snps(
  vcf_file =
    "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90017088.vcf.gz",
  target_snps = veillonella_iv$SNP,
  chunk_size = 50000
)

veillonella_eaf_raw <- veillonella_eaf_result$data
#检查
veillonella_eaf_result$n_requested
veillonella_eaf_result$n_found
veillonella_eaf_result$n_missing
veillonella_eaf_result$missing

sum(is.na(veillonella_eaf_raw$eaf))
#回填veillonellaceae_eaf_result
veillonellaceae_eaf_result <- extract_outcome_snps(
  vcf_file =
    "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90016956.vcf.gz",
  target_snps = veillonellaceae_iv$SNP,
  chunk_size = 50000
)

veillonellaceae_eaf_raw <- veillonellaceae_eaf_result$data

veillonellaceae_eaf_result$n_found

sum(is.na(veillonellaceae_eaf_raw$eaf))

# ============================================================
# 13. Identify palindromic SNPs
# ============================================================

is_palindromic <- function(a1, a2) {
  
  pair <- paste0(
    pmin(a1, a2),
    pmax(a1, a2)
  )
  
  pair %in% c("AT", "CG")
}


# Veillonella
veillonella_iv$palindromic <- is_palindromic(
  veillonella_iv$effect_allele.exposure,
  veillonella_iv$other_allele.exposure
)

veillonella_iv[, c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "palindromic"
)]

sum(veillonella_iv$palindromic)


# Veillonellaceae
veillonellaceae_iv$palindromic <- is_palindromic(
  veillonellaceae_iv$effect_allele.exposure,
  veillonellaceae_iv$other_allele.exposure
)

veillonellaceae_iv[, c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "palindromic"
)]

sum(veillonellaceae_iv$palindromic)
