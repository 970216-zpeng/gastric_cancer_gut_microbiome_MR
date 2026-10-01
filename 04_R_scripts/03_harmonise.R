# ============================================================
# 03_harmonise.R
# Project: 01_MR_GastricCancer
# Purpose:
#   Harmonise microbiome exposures with primary GC outcome
# ============================================================

source("04_R_scripts/00_functions.R")

library(TwoSampleMR)
# ============================================================
# 1. Load primary exposure IVs
# ============================================================

veillonella_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

veillonellaceae_iv <- readRDS(
  "02_intermediate_data/01_exposure/Veillonellaceae_IV_primary_p5e6.rds"
)


# ============================================================
# 2. Load primary gastric cancer outcome
# ============================================================

gc_primary <- readRDS(
  "02_intermediate_data/02_outcome/GCST90018849_primary_outcome_IVs.rds"
)
#检查读取内容
nrow(veillonella_iv)
nrow(veillonellaceae_iv)
nrow(gc_primary)

# ============================================================
# 3. Format gastric cancer outcome
# ============================================================

gc_primary$Phenotype <- "Gastric cancer"

gc_outcome <- TwoSampleMR::format_data(
  gc_primary,
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
#明确GWAS ID
gc_outcome$id.outcome <- "GCST90018849"
gc_outcome$outcome <- "Gastric cancer"
#检查
dim(gc_outcome)

head(gc_outcome)

names(gc_outcome)
#由于 exposure EAF 缺失，而我们已经决定对 palindromic SNP 采取保守策略，所以这里使用：
action = 3
# ============================================================
# 4. Harmonise Veillonella with gastric cancer
# ============================================================

harm_veillonella <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonella_iv,
  outcome_dat = gc_outcome,
  action = 3
)
#检查Harmonisation结果
nrow(harm_veillonella)
table(harm_veillonella$mr_keep)
harm_veillonella[
  ,
  c(
    "SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "palindromic",
    "ambiguous",
    "mr_keep"
  )
]
# ============================================================
# 5. Harmonise Veillonellaceae with gastric cancer
# ============================================================

harm_veillonellaceae <- TwoSampleMR::harmonise_data(
  exposure_dat = veillonellaceae_iv,
  outcome_dat = gc_outcome,
  action = 3
)
#检查
table(harm_veillonellaceae$mr_keep)
harm_veillonellaceae[
  ,
  c(
    "SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "palindromic",
    "ambiguous",
    "mr_keep"
  )
]

harm_veillonella[
  harm_veillonella$mr_keep == FALSE,
  c(
    "SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "palindromic",
    "ambiguous",
    "mr_keep"
  )
]
# ============================================================
# 6. Create MR-ready datasets
# ============================================================

veillonella_mr_ready <- subset(
  harm_veillonella,
  mr_keep == TRUE
)

veillonellaceae_mr_ready <- subset(
  harm_veillonellaceae,
  mr_keep == TRUE
)

nrow(veillonella_mr_ready)
nrow(veillonellaceae_mr_ready)
#检查
sum(is.na(veillonella_mr_ready$beta.exposure))
sum(is.na(veillonella_mr_ready$beta.outcome))
sum(is.na(veillonellaceae_mr_ready$beta.exposure))
sum(is.na(veillonellaceae_mr_ready$beta.outcome))
#保存两类文件：完整 harmonisation audit 和最终 MR-ready data 都保存
dir.create(
  "02_intermediate_data/03_harmonised",
  recursive = TRUE,
  showWarnings = FALSE
)
#完整 audit
write.csv(
  harm_veillonella,
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_harmonisation_audit.csv",
  row.names = FALSE
)

write.csv(
  harm_veillonellaceae,
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_harmonisation_audit.csv",
  row.names = FALSE
)
#完整MR-ready
write.csv(
  veillonella_mr_ready,
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_MR_ready.csv",
  row.names = FALSE
)

saveRDS(
  veillonella_mr_ready,
  "02_intermediate_data/03_harmonised/Veillonella_GCST90018849_MR_ready.rds"
)

write.csv(
  veillonellaceae_mr_ready,
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_MR_ready.csv",
  row.names = FALSE
)

saveRDS(
  veillonellaceae_mr_ready,
  "02_intermediate_data/03_harmonised/Veillonellaceae_GCST90018849_MR_ready.rds"
)
