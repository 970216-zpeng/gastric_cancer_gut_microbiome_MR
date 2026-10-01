# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})


# ============================================================
# 1. Read audited core MR master data
# ============================================================

audit_path <- file.path(
  "07_tables",
  "01_audit_data",
  "Table2_core_MR_audit_master.csv"
)

stopifnot(file.exists(audit_path))

dat <- read_csv(audit_path, show_col_types = FALSE)

stopifnot(nrow(dat) == 10)
stopifnot(all(dat$nsnp_check))

# ============================================================
# 2. Build Figure 2 plot dataset
# ============================================================

fig2_dat <- dat %>%
  mutate(
    analysis_display = case_when(
      analysis == "Primary" ~ "Primary GC",
      analysis == "Threshold sensitivity" ~ "Threshold sensitivity",
      analysis == "FinnGen sensitivity" ~ "FinnGen",
      analysis == "East Asian sensitivity" ~ "East Asian",
      analysis == "Reverse MR" ~ "Reverse MR",
      TRUE ~ analysis
    ),
    
    exposure_display = case_when(
      exposure == "Veillonella" ~ "Veillonella",
      exposure == "Veillonellaceae" ~ "Veillonellaceae",
      exposure == "Gastric cancer" ~ "Gastric cancer",
      TRUE ~ exposure
    ),
    
    panel = case_when(
      direction == "Forward" ~ "A. Forward MR",
      direction == "Reverse" ~ "B. Reverse MR"
    ),
    
    subgroup = case_when(
      direction == "Forward" & exposure == "Veillonella" ~ "Veillonella",
      direction == "Forward" & exposure == "Veillonellaceae" ~ "Veillonellaceae",
      direction == "Reverse" ~ "Reverse"
    ),
    
    y_label = case_when(
      direction == "Forward" ~ paste0(exposure_display, " — ", analysis_display),
      direction == "Reverse" ~ paste0("Gastric cancer → ", outcome)
    ),
    
    effect = case_when(
      direction == "Forward" ~ OR_recalculated_exp_beta,
      direction == "Reverse" ~ beta_output
    ),
    
    lower = case_when(
      direction == "Forward" ~ OR_CI95_low_recalculated,
      direction == "Reverse" ~ beta_ci95_low_recalculated
    ),
    
    upper = case_when(
      direction == "Forward" ~ OR_CI95_high_recalculated,
      direction == "Reverse" ~ beta_ci95_high_recalculated
    ),
    
    estimate_text = case_when(
      direction == "Forward" ~ sprintf(
        "OR %.2f (%.2f–%.2f)",
        OR_recalculated_exp_beta,
        OR_CI95_low_recalculated,
        OR_CI95_high_recalculated
      ),
      direction == "Reverse" ~ sprintf(
        "β %.4f (%.3f to %.3f)",
        beta_output,
        beta_ci95_low_recalculated,
        beta_ci95_high_recalculated
      )
    ),
    
    p_text = sprintf("%.3f", p_recalculated_from_beta_se),
    
    order_within_panel = case_when(
      y_label == "Veillonella — Primary GC" ~ 10,
      y_label == "Veillonella — Threshold sensitivity" ~ 9,
      y_label == "Veillonella — FinnGen" ~ 8,
      y_label == "Veillonella — East Asian" ~ 7,
      y_label == "Veillonellaceae — Primary GC" ~ 6,
      y_label == "Veillonellaceae — Threshold sensitivity" ~ 5,
      y_label == "Veillonellaceae — FinnGen" ~ 4,
      y_label == "Veillonellaceae — East Asian" ~ 3,
      y_label == "Gastric cancer → Veillonella" ~ 2,
      y_label == "Gastric cancer → Veillonellaceae" ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(desc(order_within_panel)) %>%
  select(
    analysis,
    analysis_display,
    direction,
    panel,
    subgroup,
    exposure,
    outcome,
    instrument_threshold,
    y_label,
    effect,
    lower,
    upper,
    estimate_text,
    p_text,
    order_within_panel
  )

# Save plotting dataset
plotdata_path <- file.path(
  "06_figures",
  "08_plot_data",
  "Figure2_core_MR_plot_data.csv"
)

write_csv(fig2_dat, plotdata_path)

stopifnot(file.exists(plotdata_path))

View(fig2_dat)
# ============================================================
# 3. Panel A: Forward MR
# ============================================================

# ============================================================
# Panel A: Forward MR
# - italicize Veillonella
# - add visual gap between the two microbial taxa
# ============================================================

fig2A_dat <- fig2_dat %>%
  filter(direction == "Forward") %>%
  mutate(
    label_md = case_when(
      exposure == "Veillonella" ~
        paste0("*Veillonella* — ", analysis_display),
      
      exposure == "Veillonellaceae" ~
        paste0("Veillonellaceae — ", analysis_display)
    ),
    
    # 手动设 y 坐标，在两个菌群之间留出间隔
    y_pos = case_when(
      exposure == "Veillonella" & analysis == "Primary" ~ 8.0,
      exposure == "Veillonella" & analysis == "Threshold sensitivity" ~ 7.0,
      exposure == "Veillonella" & analysis == "FinnGen sensitivity" ~ 6.0,
      exposure == "Veillonella" & analysis == "East Asian sensitivity" ~ 5.0,
      
      exposure == "Veillonellaceae" & analysis == "Primary" ~ 3.7,
      exposure == "Veillonellaceae" & analysis == "Threshold sensitivity" ~ 2.7,
      exposure == "Veillonellaceae" & analysis == "FinnGen sensitivity" ~ 1.7,
      exposure == "Veillonellaceae" & analysis == "East Asian sensitivity" ~ 0.7
    )
  )

pA <- ggplot(
  fig2A_dat,
  aes(x = effect, y = y_pos)
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y", 
    width = 0.16,
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.2
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = c("0.5", "0.75", "1", "1.5", "2", "3")
  ) +
  scale_y_continuous(
    breaks = fig2A_dat$y_pos,
    labels = fig2A_dat$label_md,
    limits = c(0.3, 8.5)
  ) +
  labs(
    title = "A. Forward MR",
    x = "Odds ratio for gastric cancer",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    # markdown 用于只斜体 Veillonella
    axis.text.y = ggtext::element_markdown(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
# ============================================================
# 4. Panel B: Reverse MR
# ============================================================
# ============================================================
# Panel B: Reverse MR
# - replace -> with →
# - italicize Veillonella
# ============================================================

fig2B_dat <- fig2_dat %>%
  filter(direction == "Reverse") %>%
  mutate(
    label_md = case_when(
      outcome == "Veillonella" ~
        "Gastric cancer → *Veillonella*",
      
      outcome == "Veillonellaceae" ~
        "Gastric cancer → Veillonellaceae"
    ),
    
    y_pos = case_when(
      outcome == "Veillonella" ~ 2,
      outcome == "Veillonellaceae" ~ 1
    )
  )

reverse_range <- range(
  c(fig2B_dat$lower, fig2B_dat$upper),
  na.rm = TRUE
)

reverse_pad <- 0.05 * diff(reverse_range)

pB <- ggplot(
  fig2B_dat,
  aes(x = effect, y = y_pos)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    width = 0.13,
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.2
  ) +
  scale_x_continuous(
    limits = c(
      reverse_range[1] - reverse_pad,
      reverse_range[2] + reverse_pad
    ),
    breaks = c(-0.1, 0, 0.1)
  ) +
  scale_y_continuous(
    breaks = fig2B_dat$y_pos,
    labels = fig2B_dat$label_md,
    limits = c(0.6, 2.4)
  ) +
  labs(
    title = "B. Reverse MR",
    x = expression(beta~coefficient),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text.y = ggtext::element_markdown(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
#单独保存两张图a，b
ggplot2::ggsave(
  filename = file.path("06_figures", "09_outputs", "Figure2A_forward_MR.pdf"),
  plot = pA,
  device = cairo_pdf,
  width = 7.2,
  height = 4.6,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path("06_figures", "09_outputs", "Figure2A_forward_MR.png"),
  plot = pA,
  width = 7.2,
  height = 4.6,
  dpi = 300
)

ggsave(
  filename = file.path("06_figures", "09_outputs", "Figure2B_reverse_MR.pdf"),
  plot = pB,
  device = cairo_pdf,
  width = 7.2,
  height = 2.8,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path("06_figures", "09_outputs", "Figure2B_reverse_MR.png"),
  plot = pB,
  width = 7.2,
  height = 2.8,
  dpi = 300
)

#合图调整
#合图

Figure2 <- pA / pB +
  patchwork::plot_layout(
    heights = c(4, 1.4)
  )

final_pdf <- file.path(
  "06_figures",
  "10_final",
  "Figure2_core_bidirectional_MR.pdf"
)

ggplot2::ggsave(
  filename = final_pdf,
  plot = Figure2,
  device = cairo_pdf,
  width = 7.2,
  height = 6.3,
  units = "in"
)

final_png <- file.path(
  "06_figures",
  "10_final",
  "Figure2_core_bidirectional_MR.png"
)

ggplot2::ggsave(
  filename = final_png,
  plot = Figure2,
  width = 7.2,
  height = 6.3,
  units = "in",
  dpi = 600
)

stopifnot(file.exists(final_pdf))
stopifnot(file.exists(final_png))

file.info(
  c(final_pdf, final_png)
)[, c("size", "mtime")]
#检查
file.info(c(final_pdf, final_png))[, c("size", "mtime")]

#成功标志
stopifnot(file.exists(final_pdf))
stopifnot(file.exists(final_png))

message(
  "Figure 2 script completed successfully. ",
  "Final PDF and PNG were generated."
)