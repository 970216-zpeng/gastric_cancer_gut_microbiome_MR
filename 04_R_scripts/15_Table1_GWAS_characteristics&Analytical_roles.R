# ============================================================
# Table 1. Characteristics and analytical roles of the
# GWAS datasets included in the study
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
  "Table1_publication.csv"
)

output_file <- here(
  "07_tables",
  "03_outputs",
  "Table1_GWAS_characteristics.docx"
)

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read publication data
# ------------------------------------------------------------

tab1 <- read_csv(
  input_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 3. Structural QC
# ------------------------------------------------------------

expected_cols <- c(
  "Data source",
  "Trait",
  "GWAS ID",
  "Ancestry",
  "Sample size",
  "Analytical role"
)

stopifnot(
  identical(names(tab1), expected_cols)
)

stopifnot(
  nrow(tab1) == 9
)

stopifnot(
  !anyDuplicated(tab1$`GWAS ID`)
)

# Critical metadata checks
stopifnot(
  tab1$`Sample size`[
    tab1$`GWAS ID` == "GCST90017088"
  ] == "14,306"
)

stopifnot(
  tab1$`Sample size`[
    tab1$`GWAS ID` == "GCST90016956"
  ] == "14,306"
)

stopifnot(
  tab1$`Sample size`[
    tab1$`GWAS ID` == "GCST90671339"
  ] == "16,017"
)

stopifnot(
  tab1$`Sample size`[
    tab1$`GWAS ID` == "GCST90671681"
  ] == "16,017"
)

# ------------------------------------------------------------
# 4. Build flextable
# ------------------------------------------------------------

ft <- flextable(tab1)

# General font
ft <- font(
  ft,
  fontname = "Arial",
  part = "all"
)

ft <- fontsize(
  ft,
  size = 9,
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

# Body alignment
ft <- align(
  ft,
  j = c(
    "Data source",
    "Trait",
    "GWAS ID",
    "Ancestry",
    "Analytical role"
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
# 5. Italicize genus Veillonella
# ------------------------------------------------------------

veillonella_rows <- which(
  tab1$Trait == "Veillonella"
)

for (i in veillonella_rows) {
  
  ft <- compose(
    ft,
    i = i,
    j = "Trait",
    value = as_paragraph(
      as_i("Veillonella")
    )
  )
  
}

# ------------------------------------------------------------
# 6. Booktabs-style borders
# ------------------------------------------------------------

ft <- border_remove(ft)

top_border <- fp_border(
  color = "black",
  width = 1.2
)

mid_border <- fp_border(
  color = "black",
  width = 0.7
)

bottom_border <- fp_border(
  color = "black",
  width = 1.2
)

ft <- hline_top(
  ft,
  border = top_border,
  part = "header"
)

ft <- hline_bottom(
  ft,
  border = mid_border,
  part = "header"
)

ft <- hline_bottom(
  ft,
  border = bottom_border,
  part = "body"
)

# ------------------------------------------------------------
# 7. Spacing
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

ft <- width(
  ft,
  j = "Data source",
  width = 1.65
)

ft <- width(
  ft,
  j = "Trait",
  width = 1.15
)

ft <- width(
  ft,
  j = "GWAS ID",
  width = 1.35
)

ft <- width(
  ft,
  j = "Ancestry",
  width = 0.85
)

ft <- width(
  ft,
  j = "Sample size",
  width = 1.55
)

ft <- width(
  ft,
  j = "Analytical role",
  width = 2.75
)

# Prevent individual rows breaking across pages
ft <- keep_with_next(
  ft,
  i = 1:(nrow(tab1) - 1),
  part = "body"
)

# ------------------------------------------------------------
# 9. Table title and notes
# ------------------------------------------------------------

title_text <- paste0(
  "Table 1. Characteristics and analytical roles of the ",
  "GWAS datasets"
)

note1 <- paste0(
  "Abbreviations: DMP, Dutch Microbiome Project; ",
  "GC, gastric cancer; GWAS, genome-wide association study; ",
  "MR, Mendelian randomization."
)

note2 <- paste0(
  "The MiBioGen consortium included up to 18,340 participants ",
  "across 24 cohorts; the European-ancestry summary-statistic ",
  "datasets used in the present MR analyses each comprised ",
  "14,306 participants."
)

note3 <- paste0(
  "FinnGen was treated as an alternative European outcome ",
  "sensitivity analysis rather than an independent replication ",
  "because FinnGen contributed to the primary gastric cancer ",
  "dataset GCST90018849."
)

# ------------------------------------------------------------
# 10. Build Word document
# ------------------------------------------------------------

doc <- read_docx()

# Title
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

# Table
doc <- body_add_flextable(
  doc,
  value = ft
)

# Notes
for (txt in c(note1, note2, note3)) {
  
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

# Landscape section
doc <- body_end_section_landscape(doc)

print(
  doc,
  target = output_file
)

# ------------------------------------------------------------
# 11. Final report
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("TABLE 1 GENERATED SUCCESSFULLY\n")
cat("============================================\n")
cat("Input rows: ", nrow(tab1), "\n")
cat("Input columns: ", ncol(tab1), "\n")
cat("Structural QC: PASS\n")
cat("Critical metadata QC: PASS\n")
cat("Output:\n")
cat(output_file, "\n")
cat("============================================\n")