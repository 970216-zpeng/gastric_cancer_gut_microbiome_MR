# ============================================================
# 00_functions.R
# Project: 01_MR_GastricCancer
# Purpose:
#   Reusable functions for microbiome Mendelian randomization
#
# Functions:
#   1. extract_candidate_snps()
#   2. prepare_microbiome_exposure()
#   3. extract_outcome_snps()
# ============================================================



# ============================================================
# 1. Extract candidate SNPs from microbiome GWAS-VCF
# ============================================================

extract_candidate_snps <- function(
    vcf_file,
    lp_threshold = 5,
    chunk_size = 50000
) {
  
  # Check file
  if (!file.exists(vcf_file)) {
    stop("VCF file does not exist: ", vcf_file)
  }
  
  
  con <- gzfile(
    vcf_file,
    open = "rt"
  )
  
  on.exit(
    close(con),
    add = TRUE
  )
  
  
  result_list <- list()
  result_index <- 1
  
  
  repeat {
    
    # Read VCF by chunks
    lines <- readLines(
      con,
      n = chunk_size
    )
    
    
    # End of file
    if (length(lines) == 0) {
      break
    }
    
    
    # Remove header lines
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
      
      
      # Standard GWAS-VCF should contain at least 10 columns
      if (length(fields) < 10) {
        next
      }
      
      
      # FORMAT field names
      format_names <- strsplit(
        fields[9],
        ":",
        fixed = TRUE
      )[[1]]
      
      
      # FORMAT field values
      format_values <- strsplit(
        fields[10],
        ":",
        fixed = TRUE
      )[[1]]
      
      
      # Skip malformed rows
      if (length(format_names) != length(format_values)) {
        next
      }
      
      
      names(format_values) <- format_names
      
      
      # Required fields
      if (!all(
        c("ES", "SE", "LP") %in% names(format_values)
      )) {
        next
      }
      
      
      # Convert key values
      beta <- suppressWarnings(
        as.numeric(format_values["ES"])
      )
      
      se <- suppressWarnings(
        as.numeric(format_values["SE"])
      )
      
      lp <- suppressWarnings(
        as.numeric(format_values["LP"])
      )
      
      
      if (
        is.na(beta) ||
        is.na(se) ||
        is.na(lp)
      ) {
        next
      }
      
      
      # Only retain SNPs meeting the threshold
      if (lp < lp_threshold) {
        next
      }
      
      
      # SNP ID
      if ("ID" %in% names(format_values)) {
        
        snp_id <- unname(
          format_values["ID"]
        )
        
      } else {
        
        snp_id <- fields[3]
      }
      
      
      # Alternate allele frequency
      # ES is defined relative to ALT in GWAS-VCF,
      # therefore AF corresponds to effect allele frequency.
      eaf <- if ("AF" %in% names(format_values)) {
        
        suppressWarnings(
          as.numeric(format_values["AF"])
        )
        
      } else {
        
        NA_real_
      }
      
      
      # Sample size
      samplesize <- if ("SS" %in% names(format_values)) {
        
        suppressWarnings(
          as.numeric(format_values["SS"])
        )
        
      } else {
        
        NA_real_
      }
      
      
      # Number of cases (mainly relevant for binary traits)
      ncase <- if ("NC" %in% names(format_values)) {
        
        suppressWarnings(
          as.numeric(format_values["NC"])
        )
        
      } else {
        
        NA_real_
      }
      
      
      selected_list[[selected_index]] <- data.frame(
        
        SNP = snp_id,
        
        chr = fields[1],
        
        pos = suppressWarnings(
          as.numeric(fields[2])
        ),
        
        # GWAS-VCF ES is relative to ALT allele
        effect_allele = fields[5],
        
        other_allele = fields[4],
        
        beta = beta,
        
        se = se,
        
        LP = lp,
        
        pval = 10^(-lp),
        
        eaf = eaf,
        
        samplesize = samplesize,
        
        ncase = ncase,
        
        stringsAsFactors = FALSE
      )
      
      
      selected_index <- selected_index + 1
    }
    
    
    # Save selected SNPs from this chunk
    if (length(selected_list) > 0) {
      
      result_list[[result_index]] <-
        do.call(
          rbind,
          selected_list
        )
      
      result_index <- result_index + 1
    }
  }
  
  
  # No SNPs found
  if (length(result_list) == 0) {
    
    return(
      data.frame()
    )
  }
  
  
  # Combine chunks
  result <- do.call(
    rbind,
    result_list
  )
  
  
  # Sort by P value
  result <- result[
    order(result$pval),
  ]
  
  
  # Remove duplicated rsIDs
  # Keep the strongest association
  result <- result[
    !duplicated(result$SNP),
  ]
  
  
  rownames(result) <- NULL
  
  
  return(result)
}



# ============================================================
# 2. Prepare microbiome exposure
# ============================================================

prepare_microbiome_exposure <- function(
    vcf_file,
    trait_name,
    gwas_id,
    output_prefix,
    primary_p = 5e-6,
    sensitivity_p = 1e-5,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR",
    chunk_size = 50000,
    output_dir = "02_intermediate_data/01_exposure"
) {
  
  # Check package
  if (!requireNamespace(
    "TwoSampleMR",
    quietly = TRUE
  )) {
    stop(
      "Package 'TwoSampleMR' is required."
    )
  }
  
  
  message("==========================================")
  message("Preparing exposure: ", trait_name)
  message("GWAS ID: ", gwas_id)
  message("==========================================")
  
  
  # ----------------------------------------------------------
  # Step 1. Extract SNPs using sensitivity threshold
  # ----------------------------------------------------------
  
  candidate_sensitivity <- extract_candidate_snps(
    
    vcf_file = vcf_file,
    
    lp_threshold =
      -log10(sensitivity_p),
    
    chunk_size = chunk_size
  )
  
  
  if (nrow(candidate_sensitivity) == 0) {
    
    stop(
      "No candidate SNPs found for ",
      trait_name
    )
  }
  
  
  message(
    "Candidate SNPs at P <= ",
    sensitivity_p,
    ": ",
    nrow(candidate_sensitivity)
  )
  
  
  # ----------------------------------------------------------
  # Step 2. Create primary threshold dataset
  # ----------------------------------------------------------
  
  candidate_primary <- subset(
    
    candidate_sensitivity,
    
    pval <= primary_p
  )
  
  
  if (nrow(candidate_primary) == 0) {
    
    stop(
      "No SNPs found at primary threshold for ",
      trait_name
    )
  }
  
  
  message(
    "Candidate SNPs at P <= ",
    primary_p,
    ": ",
    nrow(candidate_primary)
  )
  
  
  # ----------------------------------------------------------
  # Step 3. Calculate SNP-level F statistics
  #
  # Approximation:
  # F = (beta / SE)^2
  # ----------------------------------------------------------
  
  candidate_primary$F_stat <-
    
    (
      candidate_primary$beta /
        candidate_primary$se
    )^2
  
  
  candidate_sensitivity$F_stat <-
    
    (
      candidate_sensitivity$beta /
        candidate_sensitivity$se
    )^2
  
  
  
  # ----------------------------------------------------------
  # Step 4. Helper function for TwoSampleMR formatting
  # ----------------------------------------------------------
  
  format_exposure <- function(dat) {
    
    dat$Phenotype <- trait_name
    
    
    formatted <- TwoSampleMR::format_data(
      
      dat,
      
      type = "exposure",
      
      phenotype_col = "Phenotype",
      
      snp_col = "SNP",
      
      beta_col = "beta",
      
      se_col = "se",
      
      effect_allele_col =
        "effect_allele",
      
      other_allele_col =
        "other_allele",
      
      eaf_col = "eaf",
      
      pval_col = "pval",
      
      chr_col = "chr",
      
      pos_col = "pos",
      
      samplesize_col = "samplesize"
    )
    
    
    # Explicitly preserve GWAS ID
    formatted$id.exposure <- gwas_id
    
    formatted$exposure <- trait_name
    
    
    return(formatted)
  }
  
  
  
  # ----------------------------------------------------------
  # Step 5. Format primary exposure
  # ----------------------------------------------------------
  
  exp_primary <-
    format_exposure(
      candidate_primary
    )
  
  
  
  # ----------------------------------------------------------
  # Step 6. LD clumping: primary analysis
  # ----------------------------------------------------------
  
  clumped_primary <-
    TwoSampleMR::clump_data(
      
      exp_primary,
      
      clump_kb = clump_kb,
      
      clump_r2 = clump_r2,
      
      clump_p1 = primary_p,
      
      pop = pop
    )
  
  
  clumped_primary$F_stat <-
    
    (
      clumped_primary$beta.exposure /
        clumped_primary$se.exposure
    )^2
  
  
  
  # ----------------------------------------------------------
  # Step 7. Format sensitivity exposure
  # ----------------------------------------------------------
  
  exp_sensitivity <-
    format_exposure(
      candidate_sensitivity
    )
  
  
  
  # ----------------------------------------------------------
  # Step 8. LD clumping: threshold sensitivity analysis
  # ----------------------------------------------------------
  
  clumped_sensitivity <-
    TwoSampleMR::clump_data(
      
      exp_sensitivity,
      
      clump_kb = clump_kb,
      
      clump_r2 = clump_r2,
      
      clump_p1 = sensitivity_p,
      
      pop = pop
    )
  
  
  clumped_sensitivity$F_stat <-
    
    (
      clumped_sensitivity$beta.exposure /
        clumped_sensitivity$se.exposure
    )^2
  
  
  
  # ----------------------------------------------------------
  # Step 9. Summary table
  # ----------------------------------------------------------
  
  summary_table <- data.frame(
    
    trait = trait_name,
    
    gwas_id = gwas_id,
    
    threshold = c(
      paste0("P <= ", primary_p),
      paste0("P <= ", sensitivity_p)
    ),
    
    before_clumping = c(
      nrow(candidate_primary),
      nrow(candidate_sensitivity)
    ),
    
    after_clumping = c(
      nrow(clumped_primary),
      nrow(clumped_sensitivity)
    ),
    
    min_F = c(
      min(
        clumped_primary$F_stat,
        na.rm = TRUE
      ),
      
      min(
        clumped_sensitivity$F_stat,
        na.rm = TRUE
      )
    ),
    
    median_F = c(
      median(
        clumped_primary$F_stat,
        na.rm = TRUE
      ),
      
      median(
        clumped_sensitivity$F_stat,
        na.rm = TRUE
      )
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  
  # ----------------------------------------------------------
  # Step 10. Save results
  # ----------------------------------------------------------
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # Primary CSV
  write.csv(
    
    clumped_primary,
    
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_IV_primary_p5e6.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # Primary RDS
  saveRDS(
    
    clumped_primary,
    
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_IV_primary_p5e6.rds"
      )
    )
  )
  
  
  # Sensitivity CSV
  write.csv(
    
    clumped_sensitivity,
    
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_IV_sensitivity_p1e5.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # Sensitivity RDS
  saveRDS(
    
    clumped_sensitivity,
    
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_IV_sensitivity_p1e5.rds"
      )
    )
  )
  
  
  # Summary CSV
  write.csv(
    
    summary_table,
    
    file.path(
      output_dir,
      paste0(
        output_prefix,
        "_IV_summary.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  
  # ----------------------------------------------------------
  # Step 11. Return objects
  # ----------------------------------------------------------
  
  return(
    
    list(
      
      candidate_primary =
        candidate_primary,
      
      candidate_sensitivity =
        candidate_sensitivity,
      
      exposure_primary =
        exp_primary,
      
      exposure_sensitivity =
        exp_sensitivity,
      
      clumped_primary =
        clumped_primary,
      
      clumped_sensitivity =
        clumped_sensitivity,
      
      summary =
        summary_table
    )
  )
}



# ============================================================
# 3. Extract pre-specified SNPs from outcome GWAS-VCF
# ============================================================

extract_outcome_snps <- function(
    vcf_file,
    target_snps,
    chunk_size = 50000
) {
  
  # Check file
  if (!file.exists(vcf_file)) {
    
    stop(
      "Outcome VCF file does not exist: ",
      vcf_file
    )
  }
  
  
  target_snps <-
    unique(
      as.character(target_snps)
    )
  
  
  con <- gzfile(
    vcf_file,
    open = "rt"
  )
  
  
  on.exit(
    close(con),
    add = TRUE
  )
  
  
  result_list <- list()
  
  result_index <- 1
  
  found_snps <- character(0)
  
  
  repeat {
    
    lines <- readLines(
      con,
      n = chunk_size
    )
    
    
    if (length(lines) == 0) {
      break
    }
    
    
    # Remove VCF header
    lines <- lines[
      !startsWith(lines, "#")
    ]
    
    
    if (length(lines) == 0) {
      next
    }
    
    
    for (line in lines) {
      
      fields <- strsplit(
        line,
        "\t",
        fixed = TRUE
      )[[1]]
      
      
      if (length(fields) < 10) {
        next
      }
      
      
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
      
      
      if (
        length(format_names) !=
        length(format_values)
      ) {
        next
      }
      
      
      names(format_values) <-
        format_names
      
      
      # Obtain SNP ID
      if ("ID" %in% names(format_values)) {
        
        snp_id <-
          unname(
            format_values["ID"]
          )
        
      } else {
        
        snp_id <- fields[3]
      }
      
      
      # Skip SNPs that are not our IVs
      if (
        !(snp_id %in% target_snps)
      ) {
        next
      }
      
      
      # Required outcome values
      if (!all(
        c("ES", "SE", "LP") %in%
        names(format_values)
      )) {
        next
      }
      
      
      beta <- suppressWarnings(
        as.numeric(
          format_values["ES"]
        )
      )
      
      
      se <- suppressWarnings(
        as.numeric(
          format_values["SE"]
        )
      )
      
      
      lp <- suppressWarnings(
        as.numeric(
          format_values["LP"]
        )
      )
      
      
      # EAF
      eaf <- if (
        "AF" %in%
        names(format_values)
      ) {
        
        suppressWarnings(
          as.numeric(
            format_values["AF"]
          )
        )
        
      } else {
        
        NA_real_
      }
      
      
      # Total sample size
      samplesize <- if (
        "SS" %in%
        names(format_values)
      ) {
        
        suppressWarnings(
          as.numeric(
            format_values["SS"]
          )
        )
        
      } else {
        
        NA_real_
      }
      
      
      # Number of cases
      ncase <- if (
        "NC" %in%
        names(format_values)
      ) {
        
        suppressWarnings(
          as.numeric(
            format_values["NC"]
          )
        )
        
      } else {
        
        NA_real_
      }
      
      
      result_list[[result_index]] <-
        data.frame(
          
          SNP = snp_id,
          
          chr = fields[1],
          
          pos =
            suppressWarnings(
              as.numeric(fields[2])
            ),
          
          effect_allele =
            fields[5],
          
          other_allele =
            fields[4],
          
          beta = beta,
          
          se = se,
          
          LP = lp,
          
          pval = 10^(-lp),
          
          eaf = eaf,
          
          samplesize =
            samplesize,
          
          ncase = ncase,
          
          stringsAsFactors = FALSE
        )
      
      
      found_snps <- c(
        found_snps,
        snp_id
      )
      
      
      result_index <-
        result_index + 1
    }
    
    
    # Stop early if all target SNPs have been found
    if (
      all(
        target_snps %in%
        unique(found_snps)
      )
    ) {
      break
    }
  }
  
  
  # No SNPs found
  if (length(result_list) == 0) {
    
    result <- data.frame()
    
  } else {
    
    result <- do.call(
      rbind,
      result_list
    )
    
    
    # Keep one record per SNP
    result <- result[
      !duplicated(result$SNP),
    ]
    
    
    rownames(result) <- NULL
  }
  
  
  # Identify missing SNPs
  missing_snps <- setdiff(
    target_snps,
    result$SNP
  )
  
  
  return(
    
    list(
      
      data = result,
      
      missing = missing_snps,
      
      n_requested =
        length(target_snps),
      
      n_found =
        nrow(result),
      
      n_missing =
        length(missing_snps)
    )
  )
}


# ============================================================
# End of 00_functions.R
# ============================================================