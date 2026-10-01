############################################################
# Figure S2
# Single-SNP Mendelian randomization estimates
# Primary forward analyses
############################################################

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

############################################################
# 1. Paths
############################################################

snapshot_dir <- normalizePath(
  "../09_snapshot/2026-09-05_preregistered_core_v1.0",
  mustWork = FALSE
)

result_dir <- file.path(
  snapshot_dir,
  "05_results",
  "01_primary_MR"
)

# TODAY'S FINAL OUTPUT DIRECTORY
out_dir <- "C:/Users/zzzp1/OneDrive/Desktop/2/01_MR_GastricCancer/06_figures/10_final"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. Read single-SNP results
############################################################

v_single <- read_csv(
  file.path(
    result_dir,
    "Veillonella_GCST90018849_singleSNP.csv"
  ),
  show_col_types = FALSE
)

f_single <- read_csv(
  file.path(
    result_dir,
    "Veillonellaceae_GCST90018849_singleSNP.csv"
  ),
  show_col_types = FALSE
)


############################################################
# 3. Inspect columns
############################################################

cat("\n==============================\n")
cat("Veillonella columns\n")
cat("==============================\n")
print(names(v_single))

cat("\n==============================\n")
cat("Veillonellaceae columns\n")
cat("==============================\n")
print(names(f_single))


############################################################
# 4. Helper function
############################################################

prepare_singleSNP <- function(dat) {
  
  #------------------------------------------
  # Find SNP column
  #------------------------------------------
  
  snp_col <- intersect(
    c("SNP", "snp"),
    names(dat)
  )[1]
  
  if (is.na(snp_col)) {
    stop("SNP column not found.")
  }
  
  
  #------------------------------------------
  # Find beta and SE
  #------------------------------------------
  
  beta_col <- intersect(
    c("b", "beta"),
    names(dat)
  )[1]
  
  se_col <- intersect(
    c("se", "SE"),
    names(dat)
  )[1]
  
  if (is.na(beta_col) || is.na(se_col)) {
    stop("Cannot find beta (b) and/or SE columns.")
  }
  
  
  #------------------------------------------
  # Calculate OR and 95% CI
  #------------------------------------------
  
  dat %>%
    mutate(
      SNP_clean = as.character(.data[[snp_col]]),
      beta = as.numeric(.data[[beta_col]]),
      se = as.numeric(.data[[se_col]]),
      
      OR = exp(beta),
      
      CI_low = exp(beta - 1.96 * se),
      
      CI_high = exp(beta + 1.96 * se)
    ) %>%
    filter(
      is.finite(OR),
      is.finite(CI_low),
      is.finite(CI_high)
    )
}


############################################################
# 5. Prepare data
############################################################

v_plot <- prepare_singleSNP(v_single) %>%
  filter(!grepl("^All", SNP_clean))

f_plot <- prepare_singleSNP(f_single) %>%
  filter(!grepl("^All", SNP_clean))
############################################################
# 5A. Check raw Veillonellaceae single-SNP results
############################################################

cat("\n==============================\n")
cat("Veillonellaceae raw single-SNP results\n")
cat("==============================\n")

print(
  f_single %>%
    select(
      any_of(c("SNP", "snp")),
      any_of(c("b", "beta")),
      any_of(c("se", "SE"))
    )
)


############################################################
# 5B. Check converted Veillonellaceae values
############################################################

cat("\n==============================\n")
cat("Veillonellaceae converted values\n")
cat("==============================\n")

print(
  f_plot %>%
    select(
      SNP_clean,
      beta,
      se,
      OR,
      CI_low,
      CI_high
    )
)

############################################################
# 6. Print prepared data
############################################################

cat("\n==============================\n")
cat("Veillonella prepared data\n")
cat("==============================\n")

print(
  v_plot %>%
    select(
      SNP_clean,
      beta,
      se,
      OR,
      CI_low,
      CI_high
    )
)


cat("\n==============================\n")
cat("Veillonellaceae prepared data\n")
cat("==============================\n")

print(
  f_plot %>%
    select(
      SNP_clean,
      beta,
      se,
      OR,
      CI_low,
      CI_high
    )
)


############################################################
# 7. Forest plot function
############################################################

make_singleSNP_forest <- function(
    dat,
    panel_label,
    exposure_name
) {
  
  dat <- dat %>%
    mutate(
      SNP_clean = factor(
        SNP_clean,
        levels = rev(SNP_clean)
      )
    )
  
  
  ggplot(
    dat,
    aes(
      x = OR,
      y = SNP_clean
    )
  ) +
    
    # Reference line
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey50"
    ) +
    
    # 95% CI
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.18,
      linewidth = 0.65
    ) +
    
    # Point estimate
    geom_point(
      size = 2.3
    ) +
    
    # Log scale
    scale_x_log10() +
    
    labs(
      title = paste0(
        panel_label,
        "  ",
        exposure_name
      ),
      x = "Odds ratio (95% CI)",
      y = NULL
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      
      axis.title.x = element_text(
        size = 10.5
      ),
      
      axis.text.x = element_text(
        size = 9.5
      ),
      
      axis.text.y = element_text(
        size = 9.5
      ),
      
      axis.line = element_line(
        linewidth = 0.5
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4
      ),
      
      plot.margin = margin(
        8, 8, 8, 8
      )
    )
}


############################################################
# 8. Generate panels
############################################################

p_s2_a <- make_singleSNP_forest(
  v_plot,
  "a",
  "Veillonella"
)

p_s2_b <- make_singleSNP_forest(
  f_plot,
  "b",
  "Veillonellaceae"
)


############################################################
# 9. Combine
############################################################

fig_s2 <- (
  p_s2_a +
    p_s2_b
) +
  plot_layout(
    ncol = 2
  )


############################################################
# 10. Add global title
############################################################

fig_s2 <- fig_s2 +
  plot_annotation(
    title =
      "Single-SNP Mendelian randomization estimates for the primary forward analyses",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      )
    )
  )


############################################################
# 11. Save PDF
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Figure_S2_single_SNP_primary_forward.pdf"
  ),
  plot = fig_s2,
  width = 8.5,
  height = 4.8,
  units = "in"
)


############################################################
# 12. Save PNG
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Figure_S2_single_SNP_primary_forward.png"
  ),
  plot = fig_s2,
  width = 8.5,
  height = 4.8,
  units = "in",
  dpi = 600
)


############################################################
# 13. Confirmation
############################################################

cat("\n")
cat("==============================================\n")
cat("Figure S2 generated successfully!\n")
cat("==============================================\n")
cat(
  "PDF: ",
  file.path(
    out_dir,
    "Figure_S2_single_SNP_primary_forward.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "PNG: ",
  file.path(
    out_dir,
    "Figure_S2_single_SNP_primary_forward.png"
  ),
  "\n",
  sep = ""
)

############################################################
# Figure S3
# Leave-one-out analyses of the primary forward
# Mendelian randomization estimates
#
# DATA SOURCE:
# 05_results/01_primary_MR/*_leaveoneout.csv
#
# IMPORTANT:
# No MR analysis is recalculated in this script.
# All estimates are read directly from the frozen
# preregistered primary MR result files.
############################################################


############################################################
# 1. Libraries
############################################################

#library(readr)
#library(dplyr)
#library(ggplot2)
#library(patchwork)


############################################################
# 2. Paths
############################################################

snapshot_dir <- normalizePath(
  "../09_snapshot/2026-09-05_preregistered_core_v1.0",
  mustWork = FALSE
)

result_dir <- file.path(
  snapshot_dir,
  "05_results",
  "01_primary_MR"
)

out_dir <- "C:/Users/zzzp1/OneDrive/Desktop/2/01_MR_GastricCancer/06_figures/10_final"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Input files
############################################################

v_file <- file.path(
  result_dir,
  "Veillonella_GCST90018849_leaveoneout.csv"
)

f_file <- file.path(
  result_dir,
  "Veillonellaceae_GCST90018849_leaveoneout.csv"
)


############################################################
# 4. Check source files
############################################################

cat("\n==============================================\n")
cat("FIGURE S3 SOURCE FILE CHECK\n")
cat("==============================================\n")

cat(
  "Veillonella file exists: ",
  file.exists(v_file),
  "\n"
)

cat(
  "Veillonellaceae file exists: ",
  file.exists(f_file),
  "\n"
)

if (!file.exists(v_file)) {
  stop(
    "Veillonella leave-one-out result file not found:\n",
    v_file
  )
}

if (!file.exists(f_file)) {
  stop(
    "Veillonellaceae leave-one-out result file not found:\n",
    f_file
  )
}


############################################################
# 5. Read frozen leave-one-out results
############################################################

v_loo <- read_csv(
  v_file,
  show_col_types = FALSE
)

f_loo <- read_csv(
  f_file,
  show_col_types = FALSE
)


############################################################
# 6. Inspect source columns
############################################################

cat("\n==============================================\n")
cat("Veillonella source columns\n")
cat("==============================================\n")

print(names(v_loo))


cat("\n==============================================\n")
cat("Veillonellaceae source columns\n")
cat("==============================================\n")

print(names(f_loo))


############################################################
# 7. Helper function
############################################################

prepare_loo <- function(dat, exposure_name) {
  
  ##########################################################
  # Identify SNP column
  ##########################################################
  
  snp_col <- intersect(
    c(
      "SNP",
      "snp"
    ),
    names(dat)
  )[1]
  
  
  ##########################################################
  # Identify beta column
  ##########################################################
  
  beta_col <- intersect(
    c(
      "b",
      "beta"
    ),
    names(dat)
  )[1]
  
  
  ##########################################################
  # Identify SE column
  ##########################################################
  
  se_col <- intersect(
    c(
      "se",
      "SE"
    ),
    names(dat)
  )[1]
  
  
  if (is.na(snp_col)) {
    stop(
      exposure_name,
      ": SNP column not found."
    )
  }
  
  
  if (is.na(beta_col)) {
    stop(
      exposure_name,
      ": beta column not found."
    )
  }
  
  
  if (is.na(se_col)) {
    stop(
      exposure_name,
      ": SE column not found."
    )
  }
  
  
  ##########################################################
  # Convert beta and SE
  ##########################################################
  
  out <- dat %>%
    mutate(
      
      SNP_clean =
        as.character(.data[[snp_col]]),
      # Rename the overall estimate for publication-ready display
      SNP_clean =
        ifelse(
          SNP_clean %in% c("All", "all", "All instruments"),
          "All instruments",
          SNP_clean
        ),
      
      beta =
        as.numeric(.data[[beta_col]]),
      
      se =
        as.numeric(.data[[se_col]])
      
    )
  
  
  ##########################################################
  # Remove invalid estimates
  ##########################################################
  
  out <- out %>%
    filter(
      !is.na(SNP_clean),
      is.finite(beta),
      is.finite(se),
      se > 0
    )
  
  
  ##########################################################
  # Convert log-OR to OR
  ##########################################################
  
  out <- out %>%
    mutate(
      
      OR =
        exp(beta),
      
      CI_low =
        exp(beta - 1.96 * se),
      
      CI_high =
        exp(beta + 1.96 * se)
      
    ) %>%
    
    filter(
      is.finite(OR),
      is.finite(CI_low),
      is.finite(CI_high)
    )
  
  
  ##########################################################
  # Add exposure label
  ##########################################################
  
  out$exposure_name <- exposure_name
  
  
  return(out)
}


############################################################
# 8. Prepare datasets
############################################################

v_plot <- prepare_loo(
  v_loo,
  "Veillonella"
)

f_plot <- prepare_loo(
  f_loo,
  "Veillonellaceae"
)


############################################################
# 9. Print prepared data
############################################################

cat("\n==============================================\n")
cat("Veillonella leave-one-out estimates\n")
cat("==============================================\n")

print(
  v_plot %>%
    select(
      SNP_clean,
      beta,
      se,
      OR,
      CI_low,
      CI_high
    )
)


cat("\n==============================================\n")
cat("Veillonellaceae leave-one-out estimates\n")
cat("==============================================\n")

print(
  f_plot %>%
    select(
      SNP_clean,
      beta,
      se,
      OR,
      CI_low,
      CI_high
    )
)


############################################################
# 10. Forest plot function
############################################################

make_loo_forest <- function(
    dat,
    panel_label,
    exposure_name
) {
  
  
  ##########################################################
  # Order SNPs
  ##########################################################
  
  dat <- dat %>%
    mutate(
      SNP_clean = factor(
        SNP_clean,
        levels = rev(SNP_clean)
      )
    )
  
  
  ##########################################################
  # Plot
  ##########################################################
  
  ggplot(
    dat,
    aes(
      x = OR,
      y = SNP_clean
    )
  ) +
    
    ########################################################
  # Null effect
  ########################################################
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey50"
  ) +
    
    
    ########################################################
  # Confidence intervals
  ########################################################
  
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.18,
    linewidth = 0.65
  ) +
    
    
    ########################################################
  # Point estimates
  ########################################################
  
  geom_point(
    size = 2.3
  ) +
    
    
    ########################################################
  # Log OR scale
  ########################################################
  
  scale_x_log10() +
    
    
    ########################################################
  # Labels
  ########################################################
  
  labs(
    title =
      paste0(
        panel_label,
        "  ",
        exposure_name
      ),
    
    x =
      "Odds ratio (95% CI)",
    
    y =
      "Analysis"
  ) +
    
    
    ########################################################
  # Theme
  ########################################################
  
  theme_classic(
    base_size = 11
  ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0
        ),
      
      axis.title.x =
        element_text(
          size = 10.5
        ),
      
      axis.title.y =
        element_text(
          size = 10.5
        ),
      
      axis.text.x =
        element_text(
          size = 9.5
        ),
      
      axis.text.y =
        element_text(
          size = 9.5
        ),
      
      axis.line =
        element_line(
          linewidth = 0.5
        ),
      
      axis.ticks =
        element_line(
          linewidth = 0.4
        ),
      
      plot.margin =
        margin(
          8,
          8,
          8,
          8
        )
    )
}


############################################################
# 11. Generate panels
############################################################

p_s3_a <- make_loo_forest(
  v_plot,
  "a",
  "Veillonella"
)

p_s3_b <- make_loo_forest(
  f_plot,
  "b",
  "Veillonellaceae"
)


############################################################
# 12. Combine panels
############################################################

fig_s3 <- (
  p_s3_a +
    p_s3_b
) +
  plot_layout(
    ncol = 2
  )


############################################################
# 13. Global title
############################################################

fig_s3 <- fig_s3 +
  plot_annotation(
    title =
      "Leave-one-out analyses of the primary forward Mendelian randomization estimates",
    
    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14,
            hjust = 0
          )
      )
  )


############################################################
# 14. Save PDF
############################################################

pdf_file <- file.path(
  out_dir,
  "Figure_S3_leave_one_out_primary_forward.pdf"
)

ggsave(
  filename = pdf_file,
  plot = fig_s3,
  width = 8.5,
  height = 4.8,
  units = "in"
)


############################################################
# 15. Save PNG
############################################################

png_file <- file.path(
  out_dir,
  "Figure_S3_leave_one_out_primary_forward.png"
)

ggsave(
  filename = png_file,
  plot = fig_s3,
  width = 8.5,
  height = 4.8,
  units = "in",
  dpi = 600
)


############################################################
# 16. Export processed data for audit
############################################################

write_csv(
  v_plot,
  file.path(
    out_dir,
    "Figure_S3_Veillonella_leave_one_out_table.csv"
  )
)

write_csv(
  f_plot,
  file.path(
    out_dir,
    "Figure_S3_Veillonellaceae_leave_one_out_table.csv"
  )
)


############################################################
# 17. Audit information
############################################################

cat("\n")
cat("==============================================\n")
cat("FIGURE S3 GENERATED SUCCESSFULLY\n")
cat("==============================================\n")

cat(
  "Veillonella LOO estimates: ",
  nrow(v_plot),
  "\n"
)

cat(
  "Veillonellaceae LOO estimates: ",
  nrow(f_plot),
  "\n"
)

cat(
  "PDF: ",
  pdf_file,
  "\n",
  sep = ""
)

cat(
  "PNG: ",
  png_file,
  "\n",
  sep = ""
)

cat(
  "Processed Veillonella table: ",
  file.path(
    out_dir,
    "Figure_S3_Veillonella_leave_one_out_table.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "Processed Veillonellaceae table: ",
  file.path(
    out_dir,
    "Figure_S3_Veillonellaceae_leave_one_out_table.csv"
  ),
  "\n",
  sep = ""
)

cat("==============================================\n")

############################################################
# Figure S4
# Funnel plots of the primary forward Mendelian
# randomization analyses
############################################################

#library(readr)
#library(dplyr)
#library(ggplot2)
#library(patchwork)

############################################################
# 1. Paths
############################################################

snapshot_dir <- normalizePath(
  "../09_snapshot/2026-09-05_preregistered_core_v1.0",
  mustWork = FALSE
)

result_dir <- file.path(
  snapshot_dir,
  "05_results",
  "01_primary_MR"
)

out_dir <- "C:/Users/zzzp1/OneDrive/Desktop/2/01_MR_GastricCancer/06_figures/10_final"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. Read primary MR results
############################################################

v_mr <- read_csv(
  file.path(
    result_dir,
    "Veillonella_GCST90018849_MR.csv"
  ),
  show_col_types = FALSE
)

f_mr <- read_csv(
  file.path(
    result_dir,
    "Veillonellaceae_GCST90018849_MR.csv"
  ),
  show_col_types = FALSE
)

############################################################
# 3. Read single-SNP estimates
############################################################

v_single <- read_csv(
  file.path(
    result_dir,
    "Veillonella_GCST90018849_singleSNP.csv"
  ),
  show_col_types = FALSE
)

f_single <- read_csv(
  file.path(
    result_dir,
    "Veillonellaceae_GCST90018849_singleSNP.csv"
  ),
  show_col_types = FALSE
)

############################################################
# 4. Inspect source columns
############################################################

cat("\n==============================\n")
cat("Veillonella MR columns\n")
cat("==============================\n")
print(names(v_mr))

cat("\n==============================\n")
cat("Veillonella single-SNP columns\n")
cat("==============================\n")
print(names(v_single))

cat("\n==============================\n")
cat("Veillonellaceae MR columns\n")
cat("==============================\n")
print(names(f_mr))

cat("\n==============================\n")
cat("Veillonellaceae single-SNP columns\n")
cat("==============================\n")
print(names(f_single))

############################################################
# 5. Helper: identify columns
############################################################

find_col <- function(dat, candidates) {
  
  hit <- intersect(
    candidates,
    names(dat)
  )[1]
  
  if (is.na(hit)) {
    stop(
      "Required column not found. Candidates: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  hit
}

############################################################
# 6. Prepare single-SNP data
############################################################

prepare_funnel_data <- function(dat) {
  
  snp_col <- find_col(
    dat,
    c("SNP", "snp")
  )
  
  beta_col <- find_col(
    dat,
    c("b", "beta")
  )
  
  se_col <- find_col(
    dat,
    c("se", "SE")
  )
  
  out <- dat %>%
    mutate(
      
      SNP_clean =
        as.character(.data[[snp_col]]),
      
      beta =
        as.numeric(.data[[beta_col]]),
      
      se =
        as.numeric(.data[[se_col]])
      
    ) %>%
    
    filter(
      !is.na(SNP_clean),
      
      # Keep SNP-level estimates only
      !grepl(
        "^All\\s*-",
        SNP_clean,
        ignore.case = TRUE
      ),
      
      SNP_clean != "All",
      SNP_clean != "All instruments",
      
      is.finite(beta),
      is.finite(se),
      se > 0
    ) %>%
    
    mutate(
      precision = 1 / se
    )
  
  out
}

############################################################
# 7. Extract IVW estimate
############################################################

extract_ivw <- function(dat) {
  
  method_col <- find_col(
    dat,
    c("method", "Method")
  )
  
  beta_col <- find_col(
    dat,
    c("b", "beta")
  )
  
  ivw <- dat %>%
    filter(
      grepl(
        "Inverse variance weighted",
        .data[[method_col]],
        ignore.case = TRUE
      )
    ) %>%
    slice(1) %>%
    pull(.data[[beta_col]])
  
  if (length(ivw) == 0 || !is.finite(ivw)) {
    stop("IVW estimate not found.")
  }
  
  as.numeric(ivw)
}

############################################################
# 8. Prepare datasets
############################################################

v_funnel <- prepare_funnel_data(v_single)

f_funnel <- prepare_funnel_data(f_single)

v_ivw <- extract_ivw(v_mr)

f_ivw <- extract_ivw(f_mr)

############################################################
# 9. Print audit information
############################################################

cat("\n==============================\n")
cat("Veillonella funnel data\n")
cat("==============================\n")

print(
  v_funnel %>%
    select(
      SNP_clean,
      beta,
      se,
      precision
    )
)

cat("\nVeillonella IVW beta = ", v_ivw, "\n")

cat("\n==============================\n")
cat("Veillonellaceae funnel data\n")
cat("==============================\n")

print(
  f_funnel %>%
    select(
      SNP_clean,
      beta,
      se,
      precision
    )
)

cat("\nVeillonellaceae IVW beta = ", f_ivw, "\n")

############################################################
# 10. Funnel plot function
############################################################

make_funnel <- function(
    dat,
    ivw_beta,
    panel_label,
    exposure_name
) {
  
  ggplot(
    dat,
    aes(
      x = beta,
      y = precision
    )
  ) +
    
    # IVW reference line
    geom_vline(
      xintercept = ivw_beta,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey50"
    ) +
    
    # SNP estimates
    geom_point(
      size = 2.4
    ) +
    
    # SNP labels
    geom_text(
      aes(
        label = SNP_clean
      ),
      hjust = -0.15,
      size = 2.6
    ) +
    
    labs(
      title =
        paste0(
          panel_label,
          "  ",
          exposure_name
        ),
      
      x =
        "MR effect estimate",
      
      y =
        "1 / SE"
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0
        ),
      
      axis.title.x =
        element_text(
          size = 10.5
        ),
      
      axis.title.y =
        element_text(
          size = 10.5
        ),
      
      axis.text =
        element_text(
          size = 9.5
        ),
      
      axis.line =
        element_line(
          linewidth = 0.5
        ),
      
      axis.ticks =
        element_line(
          linewidth = 0.4
        ),
      
      plot.margin =
        margin(
          8,
          35,
          8,
          8
        )
    ) +
    
    coord_cartesian(
      clip = "off"
    )
}

############################################################
# 11. Generate panels
############################################################

p_s4_a <- make_funnel(
  v_funnel,
  v_ivw,
  "a",
  "Veillonella"
)

p_s4_b <- make_funnel(
  f_funnel,
  f_ivw,
  "b",
  "Veillonellaceae"
)

############################################################
# 12. Combine
############################################################

fig_s4 <- (
  p_s4_a +
    p_s4_b
) +
  plot_layout(
    ncol = 2
  )

############################################################
# 13. Global title
############################################################

fig_s4 <- fig_s4 +
  plot_annotation(
    title =
      "Funnel plots of the primary forward Mendelian randomization analyses",
    
    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14,
            hjust = 0
          )
      )
  )

############################################################
# 14. Save PDF
############################################################

ggsave(
  filename =
    file.path(
      out_dir,
      "Figure_S4_funnel_primary_forward.pdf"
    ),
  
  plot =
    fig_s4,
  
  width = 8.5,
  height = 4.8,
  units = "in"
)

############################################################
# 15. Save PNG
############################################################

ggsave(
  filename =
    file.path(
      out_dir,
      "Figure_S4_funnel_primary_forward.png"
    ),
  
  plot =
    fig_s4,
  
  width = 8.5,
  height = 4.8,
  units = "in",
  
  dpi = 600
)

############################################################
# 16. Save audit tables
############################################################

write_csv(
  v_funnel,
  file.path(
    out_dir,
    "Figure_S4_Veillonella_funnel_table.csv"
  )
)

write_csv(
  f_funnel,
  file.path(
    out_dir,
    "Figure_S4_Veillonellaceae_funnel_table.csv"
  )
)

############################################################
# 17. Confirmation
############################################################

cat("\n")
cat("==============================================\n")
cat("Figure S4 generated successfully!\n")
cat("==============================================\n")

cat(
  "PDF: ",
  file.path(
    out_dir,
    "Figure_S4_funnel_primary_forward.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "PNG: ",
  file.path(
    out_dir,
    "Figure_S4_funnel_primary_forward.png"
  ),
  "\n",
  sep = ""
)
############################################################
# Figure S5
# Graphical sensitivity analyses of the reverse MR analyses
#
# Row 1: Single-SNP MR estimates
# Row 2: Leave-one-out MR estimates
#
# a. Gastric cancer -> Veillonella
# b. Gastric cancer -> Veillonellaceae
# c. Gastric cancer -> Veillonella
# d. Gastric cancer -> Veillonellaceae
############################################################

#library(readr)
#library(dplyr)
#library(ggplot2)
#library(patchwork)

############################################################
# 1. Paths
############################################################

snapshot_dir <- normalizePath(
  "../09_snapshot/2026-09-05_preregistered_core_v1.0",
  mustWork = FALSE
)

result_dir <- file.path(
  snapshot_dir,
  "05_results",
  "05_reverse_MR"
)

out_dir <- "C:/Users/zzzp1/OneDrive/Desktop/2/01_MR_GastricCancer/06_figures/10_final"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 2. Read saved reverse-MR results
############################################################

v_single <- read_csv(
  file.path(
    result_dir,
    "GC_to_Veillonella_singleSNP.csv"
  ),
  show_col_types = FALSE
)

f_single <- read_csv(
  file.path(
    result_dir,
    "GC_to_Veillonellaceae_singleSNP.csv"
  ),
  show_col_types = FALSE
)

v_loo <- read_csv(
  file.path(
    result_dir,
    "GC_to_Veillonella_leaveoneout.csv"
  ),
  show_col_types = FALSE
)

f_loo <- read_csv(
  file.path(
    result_dir,
    "GC_to_Veillonellaceae_leaveoneout.csv"
  ),
  show_col_types = FALSE
)


############################################################
# 3. Inspect source columns
############################################################

cat("\n==============================\n")
cat("Reverse MR source columns\n")
cat("==============================\n\n")

cat("Veillonella Single-SNP:\n")
print(names(v_single))

cat("\nVeillonellaceae Single-SNP:\n")
print(names(f_single))

cat("\nVeillonella Leave-one-out:\n")
print(names(v_loo))

cat("\nVeillonellaceae Leave-one-out:\n")
print(names(f_loo))


############################################################
# 4. Helper function
############################################################

prepare_reverse_mr <- function(dat) {
  
  # SNP / analysis label
  snp_col <- intersect(
    c("SNP", "snp"),
    names(dat)
  )[1]
  
  if (is.na(snp_col)) {
    stop("SNP column not found.")
  }
  
  
  # beta
  beta_col <- intersect(
    c("b", "beta"),
    names(dat)
  )[1]
  
  # standard error
  se_col <- intersect(
    c("se", "SE"),
    names(dat)
  )[1]
  
  if (is.na(beta_col) || is.na(se_col)) {
    stop(
      "Cannot find beta (b/beta) and/or SE (se/SE) columns."
    )
  }
  
  
  dat %>%
    mutate(
      SNP_clean = as.character(.data[[snp_col]]),
      beta = as.numeric(.data[[beta_col]]),
      se = as.numeric(.data[[se_col]])
    ) %>%
    filter(
      is.finite(beta),
      is.finite(se)
    )
}


############################################################
# 5. Prepare data
############################################################

v_single_plot <- prepare_reverse_mr(v_single)
f_single_plot <- prepare_reverse_mr(f_single)

v_loo_plot <- prepare_reverse_mr(v_loo)
f_loo_plot <- prepare_reverse_mr(f_loo)


############################################################
# 6. Print prepared data for audit
############################################################

cat("\n============================================\n")
cat("Veillonella - Single-SNP\n")
cat("============================================\n")

print(
  v_single_plot %>%
    select(
      SNP_clean,
      beta,
      se
    )
)


cat("\n============================================\n")
cat("Veillonellaceae - Single-SNP\n")
cat("============================================\n")

print(
  f_single_plot %>%
    select(
      SNP_clean,
      beta,
      se
    )
)


cat("\n============================================\n")
cat("Veillonella - Leave-one-out\n")
cat("============================================\n")

print(
  v_loo_plot %>%
    select(
      SNP_clean,
      beta,
      se
    )
)


cat("\n============================================\n")
cat("Veillonellaceae - Leave-one-out\n")
cat("============================================\n")

print(
  f_loo_plot %>%
    select(
      SNP_clean,
      beta,
      se
    )
)

############################################################
# Plotmath title helper
# Use plotmath arrow to ensure correct PDF rendering
############################################################

#make_reverse_title <- function(
#    panel_label,
#    exposure_name
#) {
  
#  bquote(
#    bold(
#      .(panel_label) ~ "  Gastric cancer" %->% .(exposure_name)
#    )
#  )
#}
############################################################
# 7. Single-SNP forest plot
#    Only individual SNP estimates are displayed
#    Summary estimates (IVW / MR-Egger) are excluded
############################################################

make_single_plot <- function(
    dat,
    panel_label,
    exposure_name
) {
  
  #------------------------------------------
  # Keep only individual SNP estimates
  # Remove:
  #   All - Inverse variance weighted
  #   All - MR Egger
  #------------------------------------------
  
  dat <- dat %>%
    filter(
      !grepl("^All", SNP_clean, ignore.case = TRUE)
    ) %>%
    mutate(
      SNP_clean = factor(
        SNP_clean,
        levels = rev(SNP_clean)
      )
    )
  
  
  ggplot(
    dat,
    aes(
      x = beta,
      y = SNP_clean
    )
  ) +
    
    # Reference line
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey50"
    ) +
    
    # 95% CI
    geom_errorbar(
      aes(
        xmin = beta - 1.96 * se,
        xmax = beta + 1.96 * se
      ),
      orientation = "y",
      height = 0.16,
      linewidth = 0.65
    ) +
    
    # Point estimate
    geom_point(
      size = 2.2
    ) +
    
    labs(
      title = paste0(
        panel_label,
        "  Gastric cancer \u2192 ",
        exposure_name
      ),
      x = expression(
        "MR effect estimate (" * beta * ", 95% CI)"
      ),
      y = NULL
    ) +
    
    theme_classic(
      base_size = 10.5
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0
      ),
      
      axis.title.x = element_text(
        size = 10
      ),
      
      axis.text.x = element_text(
        size = 9
      ),
      
      axis.text.y = element_text(
        size = 8.5
      ),
      
      axis.line = element_line(
        linewidth = 0.5
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4
      ),
      
      plot.margin = margin(
        5, 8, 5, 5
      )
    )
}


############################################################
# 8. Leave-one-out forest plot
#    SNP labels are simplified to rsIDs
#    "Excluding" is removed
#    "All" is retained
############################################################

make_loo_plot <- function(
    dat,
    panel_label,
    exposure_name
) {
  
  dat <- dat %>%
    mutate(
      #------------------------------------------
      # Convert:
      # "Excluding rs123456"
      #      -> "rs123456"
      #
      # Keep:
      # "All"
      #------------------------------------------
      
      SNP_clean = ifelse(
        grepl(
          "^Excluding\\s+",
          SNP_clean,
          ignore.case = TRUE
        ),
        sub(
          "^Excluding\\s+",
          "",
          SNP_clean,
          ignore.case = TRUE
        ),
        SNP_clean
      )
    ) %>%
    mutate(
      SNP_clean = factor(
        SNP_clean,
        levels = rev(SNP_clean)
      )
    )
  
  
  ggplot(
    dat,
    aes(
      x = beta,
      y = SNP_clean
    )
  ) +
    
    # Reference line
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey50"
    ) +
    
    # 95% CI
    geom_errorbar(
      aes(
        xmin = beta - 1.96 * se,
        xmax = beta + 1.96 * se
      ),
      orientation = "y",
      height = 0.16,
      linewidth = 0.65
    ) +
    
    # Point estimate
    geom_point(
      size = 2.2
    ) +
    
    labs(
      title = paste0(
        panel_label,
        "  Gastric cancer \u2192 ",
        exposure_name
      ),
      x = expression(
        "Leave-one-out MR estimate (" * beta * ", 95% CI)"
      ),
      y = NULL
    ) +
    
    theme_classic(
      base_size = 10.5
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        hjust = 0
      ),
      
      axis.title.x = element_text(
        size = 10
      ),
      
      axis.text.x = element_text(
        size = 9
      ),
      
      axis.text.y = element_text(
        size = 8.5
      ),
      
      axis.line = element_line(
        linewidth = 0.5
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4
      ),
      
      plot.margin = margin(
        5, 8, 5, 5
      )
    )
}

############################################################
# 9. Generate four panels
############################################################

p_s5_a <- make_single_plot(
  v_single_plot,
  "a",
  "Veillonella"
)

p_s5_b <- make_single_plot(
  f_single_plot,
  "b",
  "Veillonellaceae"
)

p_s5_c <- make_loo_plot(
  v_loo_plot,
  "c",
  "Veillonella"
)

p_s5_d <- make_loo_plot(
  f_loo_plot,
  "d",
  "Veillonellaceae"
)


############################################################
# 10. Combine 2 × 2
############################################################

fig_s5 <- (
  (p_s5_a | p_s5_b) /
    (p_s5_c | p_s5_d)
)


############################################################
# 11. Global title
############################################################

fig_s5 <- fig_s5 +
  plot_annotation(
    title =
      "Graphical sensitivity analyses of the reverse Mendelian randomization analyses",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      )
    )
  )


############################################################
# 12. Save PDF
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Figure_S5_reverse_MR_sensitivity.pdf"
  ),
  plot = fig_s5,
  width = 8.5,
  height = 7.2,
  units = "in",
  device = cairo_pdf
)


############################################################
# 13. Save PNG
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Figure_S5_reverse_MR_sensitivity.png"
  ),
  plot = fig_s5,
  width = 8.5,
  height = 7.2,
  units = "in",
  dpi = 600
)


############################################################
# 14. Confirmation
############################################################

cat("\n")
cat("==============================================\n")
cat("Figure S5 generated successfully!\n")
cat("==============================================\n")

cat(
  "PDF: ",
  file.path(
    out_dir,
    "Figure_S5_reverse_MR_sensitivity.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "PNG: ",
  file.path(
    out_dir,
    "Figure_S5_reverse_MR_sensitivity.png"
  ),
  "\n",
  sep = ""
)
############################################################
# Figure S6
# Instrument-level sensitivity analyses of the
# secondary-threshold DMP Veillonellaceae association
# with gastric cancer
############################################################


############################################################
# 1. Packages
############################################################

library(tidyverse)
library(patchwork)


############################################################
# 2. File paths
############################################################

data_dir <- paste0(
  "C:/Users/zzzp1/OneDrive/Desktop/2/",
  "01_MR_GastricCancer/09_snapshot/",
  "2026-09-07_postregistration_exposure_replication_v1.0/",
  "05_results/07_exposure_source_replication"
)

single_file <- file.path(
  data_dir,
  "external_singleSNP.csv"
)

loo_file <- file.path(
  data_dir,
  "external_leaveoneout.csv"
)

out_dir <- data_dir


############################################################
# 3. Read raw data
############################################################

external_single <- read.csv(
  single_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

external_loo <- read.csv(
  loo_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 4. Target definition
############################################################

target_exposure_id <- "GCST90027679"

target_analysis <-
  "DMP2022_Veillonellaceae_sensitivity_p1e5"


############################################################
# 5. RAW DATA AUDIT
############################################################

cat("\n")
cat("====================================================\n")
cat("S6 RAW DATA AUDIT\n")
cat("====================================================\n")

cat(
  "\nSingle-SNP raw rows: ",
  nrow(external_single),
  "\n",
  sep = ""
)

cat(
  "LOO raw rows: ",
  nrow(external_loo),
  "\n",
  sep = ""
)


############################################################
# 6. Extract secondary-threshold Single-SNP data
############################################################

s6_single <- external_single %>%
  
  filter(
    id.exposure == target_exposure_id,
    analysis == target_analysis
  ) %>%
  
  mutate(
    SNP = trimws(as.character(SNP)),
    b = as.numeric(b),
    se = as.numeric(se)
  ) %>%
  
  # Keep actual SNP-level estimates only
  filter(
    !grepl("^All", SNP, ignore.case = TRUE)
  ) %>%
  
  mutate(
    OR = exp(b),
    
    CI_low = exp(
      b - 1.96 * se
    ),
    
    CI_high = exp(
      b + 1.96 * se
    )
  ) %>%
  
  filter(
    is.finite(OR),
    is.finite(CI_low),
    is.finite(CI_high)
  )


############################################################
# 7. Extract secondary-threshold LOO data
############################################################

s6_loo <- external_loo %>%
  
  filter(
    id.exposure == target_exposure_id,
    analysis == target_analysis
  ) %>%
  
  mutate(
    SNP = trimws(as.character(SNP)),
    b = as.numeric(b),
    se = as.numeric(se)
  ) %>%
  
  mutate(
    SNP = ifelse(
      grepl(
        "^All",
        SNP,
        ignore.case = TRUE
      ),
      "All instruments",
      SNP
    ),
    
    OR = exp(b),
    
    CI_low = exp(
      b - 1.96 * se
    ),
    
    CI_high = exp(
      b + 1.96 * se
    )
  ) %>%
  
  filter(
    is.finite(OR),
    is.finite(CI_low),
    is.finite(CI_high)
  )


############################################################
# 8. DATA AUDIT
############################################################

cat("\n")
cat("====================================================\n")
cat("S6 FILTERED DATA AUDIT\n")
cat("====================================================\n")

cat(
  "\nSingle-SNP rows: ",
  nrow(s6_single),
  "\n",
  sep = ""
)

cat(
  "LOO rows: ",
  nrow(s6_loo),
  "\n",
  sep = ""
)


cat("\nSingle-SNP SNPs:\n")
print(s6_single$SNP)


cat("\nLOO SNPs / labels:\n")
print(s6_loo$SNP)


cat("\nSingle-SNP OR range:\n")
print(range(s6_single$OR))


cat("\nLOO OR range:\n")
print(range(s6_loo$OR))


############################################################
# 9. Stop if extraction failed
############################################################

if (nrow(s6_single) == 0) {
  
  stop(
    "\nERROR: No Single-SNP records were extracted."
  )
}


if (nrow(s6_loo) == 0) {
  
  stop(
    "\nERROR: No leave-one-out records were extracted."
  )
}


############################################################
# 10. Prepare plotting order
############################################################

s6_single <- s6_single %>%
  
  mutate(
    SNP = factor(
      SNP,
      levels = rev(SNP)
    )
  )


s6_loo <- s6_loo %>%
  
  mutate(
    SNP = factor(
      SNP,
      levels = rev(
        unique(SNP)
      )
    )
  )


############################################################
# 11. Single-SNP forest plot
############################################################

p_s6_a <- ggplot(
  s6_single,
  aes(
    x = OR,
    y = SNP
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey50"
  ) +
  
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    orientation = "y",
    height = 0.14,
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 2.3
  ) +
  
  labs(
    title = "a  Single-SNP estimates",
    x = "Odds ratio (95% CI)",
    y = NULL
  ) +
  
  scale_x_continuous(
    limits = c(
      min(s6_single$CI_low) * 0.95,
      max(s6_single$CI_high) * 1.05
    ),
    expand = expansion(
      mult = c(0.02, 0.03)
    )
  ) +
  
  theme_classic(
    base_size = 10.5
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 11.5,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 8.5
    ),
    
    axis.line = element_line(
      linewidth = 0.5
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      5, 8, 5, 5
    )
  )


############################################################
# 12. Leave-one-out forest plot
############################################################

p_s6_b <- ggplot(
  s6_loo,
  aes(
    x = OR,
    y = SNP
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey50"
  ) +
  
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    orientation = "y",
    height = 0.14,
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 2.3
  ) +
  
  labs(
    title = "b  Leave-one-out estimates",
    x = "Odds ratio (95% CI)",
    y = NULL
  ) +
  
  scale_x_continuous(
    limits = c(
      min(s6_loo$CI_low) * 0.95,
      max(s6_loo$CI_high) * 1.05
    ),
    expand = expansion(
      mult = c(0.02, 0.03)
    )
  ) +
  
  theme_classic(
    base_size = 10.5
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 11.5,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 8.5
    ),
    
    axis.line = element_line(
      linewidth = 0.5
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      5, 8, 5, 5
    )
  )


############################################################
# 13. Combine panels
############################################################

fig_s6 <-
  p_s6_a | p_s6_b


############################################################
# 14. Global title
############################################################

#fig_s6 <-
#  fig_s6 +
  
#  plot_annotation(
    
#    title = "null",
    
#    theme =
#      theme(
        
#        plot.title =
#          element_text(
#            face = "bold",
#            size = 14,
#            hjust = 0
#          )
#      )
#  )


############################################################
# 15. Save PDF
############################################################

ggsave(
  
  filename =
    file.path(
      out_dir,
      "Figure_S6_DMP_Veillonellaceae_secondary_threshold_sensitivity.pdf"
    ),
  
  plot = fig_s6,
  
  width = 10,
  height = 5.2,
  units = "in"
)


############################################################
# 16. Save PNG
############################################################

ggsave(
  
  filename =
    file.path(
      out_dir,
      "Figure_S6_DMP_Veillonellaceae_secondary_threshold_sensitivity.png"
    ),
  
  plot = fig_s6,
  
  width = 10,
  height = 5.2,
  units = "in",
  
  dpi = 600
)


############################################################
# 17. Final confirmation
############################################################

cat("\n")
cat("====================================================\n")
cat("Figure S6 generated successfully!\n")
cat("====================================================\n")

cat(
  "\nPDF:\n",
  file.path(
    out_dir,
    "Figure_S6_DMP_Veillonellaceae_secondary_threshold_sensitivity.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "\nPNG:\n",
  file.path(
    out_dir,
    "Figure_S6_DMP_Veillonellaceae_secondary_threshold_sensitivity.png"
  ),
  "\n",
  sep = ""
)