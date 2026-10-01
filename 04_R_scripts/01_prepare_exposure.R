# ============================================================
# 01_prepare_exposure.R
# Project: 01_MR_GastricCancer
# Purpose:
#   Prepare MiBioGen exposure GWAS
#   Primary taxon: Veillonella
# ============================================================


# ============================================================
# 0. Check working directory
# ============================================================

getwd()

# 应该类似：
# C:/Users/zzzp1/OneDrive/Desktop/2/01_MR_GastricCancer


# ============================================================
# 1. Define file paths
# ============================================================

# Family: Veillonellaceae
file_veillonellaceae <-
  "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90016956.vcf.gz"

# Genus: Veillonella
file_veillonella <-
  "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90017088.vcf.gz"


# ============================================================
# 2. Check files
# ============================================================

file.exists(file_veillonellaceae)
file.exists(file_veillonella)

# 两个都应该返回 TRUE


# 查看文件大小（MB）
file.info(file_veillonellaceae)$size / 1024^2
file.info(file_veillonella)$size / 1024^2
# ============================================================
# 3. Read Veillonella VCF header
# ============================================================

con <- gzfile(file_veillonella, open = "rt")

vcf_head <- readLines(con, n = 200)

close(con)


# 查看前 20 行
vcf_head[1:20]


# 找 VCF 列名
vcf_head[grepl("^#CHROM", vcf_head)]


# 查看 FORMAT 定义
vcf_head[grepl("^##FORMAT", vcf_head)]
# ============================================================
# 4. Read first SNP record，确认解析正常
# ============================================================

con <- gzfile(file_veillonella, open = "rt")

repeat {
  
  line <- readLines(con, n = 1)
  
  if (length(line) == 0) {
    break
  }
  
  if (!startsWith(line, "#")) {
    
    first_snp_line <- line
    break
  }
}

close(con)


# 查看原始 SNP 行
first_snp_line
# ============================================================
# 5. Parse first SNP record
# ============================================================

# 获取列名
header_line <- vcf_head[grepl("^#CHROM", vcf_head)]

vcf_colnames <- strsplit(
  sub("^#", "", header_line),
  "\t"
)[[1]]

vcf_colnames


# 拆分 SNP 行
snp_fields <- strsplit(
  first_snp_line,
  "\t",
  fixed = TRUE
)[[1]]

names(snp_fields) <- vcf_colnames


# 看 FORMAT
snp_fields["FORMAT"]


# 拆 FORMAT 名称
format_names <- strsplit(
  snp_fields["FORMAT"],
  ":",
  fixed = TRUE
)[[1]]


# 拆 FORMAT 数值
format_values <- strsplit(
  snp_fields["ebi-a-GCST90017088"],
  ":",
  fixed = TRUE
)[[1]]

names(format_values) <- format_names


format_values
# ============================================================
# 6. Convert first SNP into MR-style format
# ============================================================

first_snp <- data.frame(
  
  SNP = unname(format_values["ID"]),
  
  chr = unname(snp_fields["CHROM"]),
  
  pos = as.numeric(snp_fields["POS"]),
  
  effect_allele = unname(snp_fields["ALT"]),
  
  other_allele = unname(snp_fields["REF"]),
  
  beta = as.numeric(format_values["ES"]),
  
  se = as.numeric(format_values["SE"]),
  
  LP = as.numeric(format_values["LP"]),
  
  pval = 10^(-as.numeric(format_values["LP"])),
  
  stringsAsFactors = FALSE
)

first_snp
# ============================================================
# 7. Function: extract candidate SNPs from GWAS-VCF
# ============================================================

extract_candidate_snps <- function(
    vcf_file,
    lp_threshold = 5,
    chunk_size = 50000
) {
  
  con <- gzfile(vcf_file, open = "rt")
  
  on.exit(close(con))
  
  result_list <- list()
  
  result_index <- 1
  
  repeat {
    
    # 每次读取一部分，避免一次把整个 VCF 放进内存
    lines <- readLines(
      con,
      n = chunk_size
    )
    
    # 文件读完
    if (length(lines) == 0) {
      break
    }
    
    # 去掉 VCF header
    lines <- lines[
      !startsWith(lines, "#")
    ]
    
    if (length(lines) == 0) {
      next
    }
    
    
    selected_list <- list()
    
    selected_index <- 1
    
    
    for (line in lines) {
      
      fields <- strsplit(
        line,
        "\t",
        fixed = TRUE
      )[[1]]
      
      
      # 至少应该有 10 列
      if (length(fields) < 10) {
        next
      }
      
      
      # FORMAT，例如：
      # ES:SE:LP:ID
      format_names <- strsplit(
        fields[9],
        ":",
        fixed = TRUE
      )[[1]]
      
      
      format_values <- strsplit(
        fields[10],
        ":",
        fixed = TRUE
      )[[1]]
      
      
      # 名称和值数量不一致则跳过
      if (length(format_names) != length(format_values)) {
        next
      }
      
      
      names(format_values) <- format_names
      
      
      # 必须至少有 ES、SE、LP
      if (!all(c("ES", "SE", "LP") %in% names(format_values))) {
        next
      }
      
      
      lp <- suppressWarnings(
        as.numeric(format_values["LP"])
      )
      
      
      # LP 缺失则跳过
      if (is.na(lp)) {
        next
      }
      
      
      # ==========================================
      # Candidate IV selection:
      # LP >= 5 corresponds to P <= 1e-5
      # ==========================================
      
      if (lp >= lp_threshold) {
        
        
        # 优先使用 FORMAT 中的 rsID
        if ("ID" %in% names(format_values)) {
          
          snp_id <- format_values["ID"]
          
        } else {
          
          snp_id <- fields[3]
          
        }
        
        
        selected_list[[selected_index]] <- data.frame(
          
          SNP = unname(snp_id),
          
          chr = fields[1],
          
          pos = as.numeric(fields[2]),
          
          effect_allele = fields[5],
          
          other_allele = fields[4],
          
          beta = as.numeric(format_values["ES"]),
          
          se = as.numeric(format_values["SE"]),
          
          LP = lp,
          
          pval = 10^(-lp),
          
          stringsAsFactors = FALSE
        )
        
        
        selected_index <- selected_index + 1
      }
    }
    
    
    # 保存这一批中筛到的 SNP
    if (length(selected_list) > 0) {
      
      result_list[[result_index]] <-
        do.call(
          rbind,
          selected_list
        )
      
      result_index <- result_index + 1
    }
  }
  
  
  # 如果一个 SNP 都没找到
  if (length(result_list) == 0) {
    
    return(data.frame())
  }
  
  
  # 合并所有批次
  result <- do.call(
    rbind,
    result_list
  )
  
  
  rownames(result) <- NULL
  
  return(result)
}
# ============================================================
# 8. Extract Veillonella candidate instruments
#    Threshold: P <= 1e-5
# ============================================================

veillonella_p1e5 <- extract_candidate_snps(
  
  vcf_file = file_veillonella,
  
  lp_threshold = 5,
  
  chunk_size = 50000
)
# ============================================================
# 9. Inspect candidate instruments
# ============================================================

dim(veillonella_p1e5)

nrow(veillonella_p1e5)

head(veillonella_p1e5)

summary(veillonella_p1e5$pval)
#检查最大和最小Pval
min(veillonella_p1e5$pval, na.rm = TRUE)

max(veillonella_p1e5$pval, na.rm = TRUE)
# ============================================================
# 10. Compare different P-value thresholds
# ============================================================

# P <= 1e-5
veillonella_p1e5_final <- veillonella_p1e5

# P <= 5e-6
veillonella_p5e6 <- subset(
  veillonella_p1e5,
  pval <= 5e-6
)

# Compare number of candidate SNPs
nrow(veillonella_p1e5_final)
nrow(veillonella_p5e6)
# ============================================================
# 11. Check chromosome distribution
# ============================================================

# P <= 1e-5
table(veillonella_p1e5$chr)

# P <= 5e-6
table(veillonella_p5e6$chr)
#检查严格阈值下16个SNP的位置
veillonella_p5e6[
  order(
    as.numeric(veillonella_p5e6$chr),
    veillonella_p5e6$pos
  ),
  c(
    "SNP",
    "chr",
    "pos",
    "beta",
    "se",
    "pval"
  )
]
# ============================================================
# 12. Calculate F statistics
# ============================================================

veillonella_p1e5$F_stat <-
  (veillonella_p1e5$beta / veillonella_p1e5$se)^2

veillonella_p5e6$F_stat <-
  (veillonella_p5e6$beta / veillonella_p5e6$se)^2
#先看看P1e5
summary(veillonella_p1e5$F_stat)

min(veillonella_p1e5$F_stat)

sum(veillonella_p1e5$F_stat <= 10)
#再看P5e6
summary(veillonella_p5e6$F_stat)

min(veillonella_p5e6$F_stat)

sum(veillonella_p5e6$F_stat <= 10)
#最后可以查看最弱的几个，如果有
veillonella_p1e5[
  order(veillonella_p1e5$F_stat),
  c("SNP", "beta", "se", "pval", "F_stat")
][1:10, ]

# ============================================================
# 13. Load TwoSampleMR
# ============================================================
install.packages(
  "TwoSampleMR",
  repos = c(
    "https://mrcieu.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
library(TwoSampleMR)

packageVersion("TwoSampleMR")

# ============================================================
# 14. Format exposure data for TwoSampleMR
# ============================================================

veillonella_p5e6$Phenotype <- "Veillonella"

veillonella_exp_p5e6 <- format_data(
  
  veillonella_p5e6,
  
  type = "exposure",
  
  phenotype_col = "Phenotype",
  
  snp_col = "SNP",
  
  beta_col = "beta",
  
  se_col = "se",
  
  effect_allele_col = "effect_allele",
  
  other_allele_col = "other_allele",
  
  pval_col = "pval",
  
  chr_col = "chr",
  
  pos_col = "pos"
)
#检查
dim(veillonella_exp_p5e6)

head(veillonella_exp_p5e6)

names(veillonella_exp_p5e6)

# ============================================================
# 15. LD clumping
# ============================================================

veillonella_clumped_p5e6 <- clump_data(
  
  veillonella_exp_p5e6,
  
  clump_kb = 10000,
  
  clump_r2 = 0.001,
  
  clump_p1 = 5e-6,
  
  pop = "EUR"
)
# ============================================================
#OpenGWAS API 需要 JWT token，去创建账号，创token

#install.packages("usethis")   # 如果以前装过可跳过
#library(usethis)

#usethis::edit_r_environ()
#重启R，session-重启；并验证R是否读到token

#library(ieugwasr)

#jwt <- ieugwasr::get_opengwas_jwt()

#nchar(jwt)
#再测试认证，返回正常账户，则JWT配置成功
#ieugwasr::user()
# ============================================================

# ============================================================
# 16. Inspect clumped Veillonella instruments
# ============================================================

nrow(veillonella_clumped_p5e6)

veillonella_clumped_p5e6[
  ,
  c(
    "SNP",
    "chr.exposure",
    "pos.exposure",
    "beta.exposure",
    "se.exposure",
    "pval.exposure"
  )
]
# Add F-statistic back to clumped IVs

veillonella_clumped_p5e6$F_stat <-
  (veillonella_clumped_p5e6$beta.exposure /
     veillonella_clumped_p5e6$se.exposure)^2

veillonella_clumped_p5e6[
  ,
  c(
    "SNP",
    "chr.exposure",
    "pos.exposure",
    "beta.exposure",
    "se.exposure",
    "pval.exposure",
    "F_stat"
  )
]
#再检查
summary(veillonella_clumped_p5e6$F_stat)

min(veillonella_clumped_p5e6$F_stat)

# ============================================================
# 17. Format Veillonella exposure: P <= 1e-5
# ============================================================

veillonella_p1e5$Phenotype <- "Veillonella"

veillonella_exp_p1e5 <- format_data(
  veillonella_p1e5,
  type = "exposure",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "pos"
)

dim(veillonella_exp_p1e5)
# ============================================================
# 18. LD clumping: P <= 1e-5
# ============================================================

veillonella_clumped_p1e5 <- clump_data(
  veillonella_exp_p1e5,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1e-5,
  pop = "EUR"
)

nrow(veillonella_clumped_p1e5)

veillonella_clumped_p1e5$SNP
#补上F值
veillonella_clumped_p1e5$F_stat <-
  (veillonella_clumped_p1e5$beta.exposure /
     veillonella_clumped_p1e5$se.exposure)^2

summary(veillonella_clumped_p1e5$F_stat)

min(veillonella_clumped_p1e5$F_stat)

#比较两个阈值
# ============================================================
# 19. Compare IV sets
# ============================================================

data.frame(
  threshold = c("P <= 5e-6", "P <= 1e-5"),
  before_clumping = c(
    nrow(veillonella_p5e6),
    nrow(veillonella_p1e5)
  ),
  after_clumping = c(
    nrow(veillonella_clumped_p5e6),
    nrow(veillonella_clumped_p1e5)
  ),
  min_F = c(
    min(veillonella_clumped_p5e6$F_stat),
    min(veillonella_clumped_p1e5$F_stat)
  )
)
#查看阈值宽松后，增加的SNP
setdiff(
  veillonella_clumped_p1e5$SNP,
  veillonella_clumped_p5e6$SNP
)
# ============================================================
# 20. Save Veillonella IV results
# ============================================================

dir.create(
  "02_intermediate_data/01_exposure",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  veillonella_clumped_p5e6,
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.csv",
  row.names = FALSE
)

write.csv(
  veillonella_clumped_p1e5,
  "02_intermediate_data/01_exposure/Veillonella_IV_sensitivity_p1e5.csv",
  row.names = FALSE
)

saveRDS(
  veillonella_clumped_p5e6,
  "02_intermediate_data/01_exposure/Veillonella_IV_primary_p5e6.rds"
)

saveRDS(
  veillonella_clumped_p1e5,
  "02_intermediate_data/01_exposure/Veillonella_IV_sensitivity_p1e5.rds"
)

# ============================================================
# 21. Reusable function for microbiome exposure preparation
# ============================================================

prepare_microbiome_exposure <- function(
    vcf_file,
    trait_name,
    gwas_id,
    output_prefix
) {
  
  message("==========================================")
  message("Preparing exposure: ", trait_name)
  message("GWAS ID: ", gwas_id)
  message("==========================================")
  
  
  # ----------------------------------------------------------
  # Step 1. Extract candidate SNPs at P <= 1e-5
  # ----------------------------------------------------------
  
  candidate_p1e5 <- extract_candidate_snps(
    vcf_file = vcf_file,
    lp_threshold = 5,
    chunk_size = 50000
  )
  
  message(
    "Candidate SNPs at P <= 1e-5: ",
    nrow(candidate_p1e5)
  )
  
  
  # ----------------------------------------------------------
  # Step 2. Create stricter P <= 5e-6 dataset
  # ----------------------------------------------------------
  
  candidate_p5e6 <- subset(
    candidate_p1e5,
    pval <= 5e-6
  )
  
  message(
    "Candidate SNPs at P <= 5e-6: ",
    nrow(candidate_p5e6)
  )
  
  
  # ----------------------------------------------------------
  # Step 3. Calculate F statistics before clumping
  # ----------------------------------------------------------
  
  candidate_p1e5$F_stat <-
    (candidate_p1e5$beta / candidate_p1e5$se)^2
  
  candidate_p5e6$F_stat <-
    (candidate_p5e6$beta / candidate_p5e6$se)^2
  
  
  # ----------------------------------------------------------
  # Step 4. Format P <= 5e-6 data for TwoSampleMR
  # ----------------------------------------------------------
  
  candidate_p5e6$Phenotype <- trait_name
  
  exp_p5e6 <- format_data(
    candidate_p5e6,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )
  
  
  # ----------------------------------------------------------
  # Step 5. LD clumping: primary threshold P <= 5e-6
  # ----------------------------------------------------------
  
  clumped_p5e6 <- clump_data(
    exp_p5e6,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 5e-6,
    pop = "EUR"
  )
  
  clumped_p5e6$F_stat <-
    (clumped_p5e6$beta.exposure /
       clumped_p5e6$se.exposure)^2
  
  
  # ----------------------------------------------------------
  # Step 6. Format P <= 1e-5 data
  # ----------------------------------------------------------
  
  candidate_p1e5$Phenotype <- trait_name
  
  exp_p1e5 <- format_data(
    candidate_p1e5,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )
  
  
  # ----------------------------------------------------------
  # Step 7. LD clumping: sensitivity threshold P <= 1e-5
  # ----------------------------------------------------------
  
  clumped_p1e5 <- clump_data(
    exp_p1e5,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 1e-5,
    pop = "EUR"
  )
  
  clumped_p1e5$F_stat <-
    (clumped_p1e5$beta.exposure /
       clumped_p1e5$se.exposure)^2
  
  
  # ----------------------------------------------------------
  # Step 8. Summary table
  # ----------------------------------------------------------
  
  summary_table <- data.frame(
    
    trait = trait_name,
    
    threshold = c(
      "P <= 5e-6",
      "P <= 1e-5"
    ),
    
    before_clumping = c(
      nrow(candidate_p5e6),
      nrow(candidate_p1e5)
    ),
    
    after_clumping = c(
      nrow(clumped_p5e6),
      nrow(clumped_p1e5)
    ),
    
    min_F = c(
      min(clumped_p5e6$F_stat),
      min(clumped_p1e5$F_stat)
    )
  )
  
  
  # ----------------------------------------------------------
  # Step 9. Save results
  # ----------------------------------------------------------
  
  dir.create(
    "02_intermediate_data/01_exposure",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  write.csv(
    clumped_p5e6,
    paste0(
      "02_intermediate_data/01_exposure/",
      output_prefix,
      "_IV_primary_p5e6.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    clumped_p1e5,
    paste0(
      "02_intermediate_data/01_exposure/",
      output_prefix,
      "_IV_sensitivity_p1e5.csv"
    ),
    row.names = FALSE
  )
  
  
  saveRDS(
    clumped_p5e6,
    paste0(
      "02_intermediate_data/01_exposure/",
      output_prefix,
      "_IV_primary_p5e6.rds"
    )
  )
  
  
  saveRDS(
    clumped_p1e5,
    paste0(
      "02_intermediate_data/01_exposure/",
      output_prefix,
      "_IV_sensitivity_p1e5.rds"
    )
  )
  
  
  # ----------------------------------------------------------
  # Step 10. Return all useful objects
  # ----------------------------------------------------------
  
  return(
    list(
      
      candidate_p5e6 = candidate_p5e6,
      
      candidate_p1e5 = candidate_p1e5,
      
      clumped_p5e6 = clumped_p5e6,
      
      clumped_p1e5 = clumped_p1e5,
      
      summary = summary_table
    )
  )
}
# ============================================================
# 22. Prepare Veillonellaceae
# ============================================================

veillonellaceae_result <- prepare_microbiome_exposure(
  
  vcf_file =
    "01_raw_data/01_exposure/MiBioGen/ebi-a-GCST90016956.vcf.gz",
  
  trait_name =
    "Veillonellaceae",
  
  gwas_id =
    "GCST90016956",
  
  output_prefix =
    "Veillonellaceae"
)
#检查
veillonellaceae_result$summary

# ============================================================
# 23. Compare genus-level and family-level IVs
# ============================================================

intersect(
  veillonella_clumped_p5e6$SNP,
  veillonellaceae_result$clumped_p5e6$SNP
)
length(
  intersect(
    veillonella_clumped_p5e6$SNP,
    veillonellaceae_result$clumped_p5e6$SNP
  )
)
# Veillonella-specific IVs
setdiff(
  veillonella_clumped_p5e6$SNP,
  veillonellaceae_result$clumped_p5e6$SNP
)

# Veillonellaceae-specific IVs
setdiff(
  veillonellaceae_result$clumped_p5e6$SNP,
  veillonella_clumped_p5e6$SNP
)
