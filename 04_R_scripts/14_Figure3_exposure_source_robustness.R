# ============================================================
# Figure 3: Exposure-source robustness
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

# ============================================================
# 1. Read independently audited Figure 3 master data
# ============================================================

audit_path <- file.path(
  "07_tables",
  "01_audit_data",
  "Figure3_exposure_source_audit_master.csv"
)

stopifnot(file.exists(audit_path))

dat <- readr::read_csv(
  audit_path,
  show_col_types = FALSE
)

# Basic QC
stopifnot(nrow(dat) == 12)
stopifnot(all(dat$nsnp_check))
# ============================================================
# 2. Build Figure 3 plotting dataset
# ============================================================

fig3_dat <- dat %>%
  mutate(
    
    threshold_display = case_when(
      instrument_threshold == "P <= 5e-6" ~
        "P ≤ 5 × 10⁻⁶",
      
      instrument_threshold == "P <= 1e-5" ~
        "P ≤ 1 × 10⁻⁵",
      
      TRUE ~ instrument_threshold
    ),
    
    threshold_role = case_when(
      analysis_role == "Primary threshold" ~ "Primary",
      analysis_role == "Secondary threshold sensitivity" ~ "Secondary",
      TRUE ~ analysis_role
    ),
    
    source_order = case_when(
      source == "MiBioGen" ~ 1,
      source == "Swedish 2026" ~ 2,
      source == "DMP 2022" ~ 3
    ),
    
    threshold_order = case_when(
      analysis_role == "Primary threshold" ~ 1,
      analysis_role == "Secondary threshold sensitivity" ~ 2
    ),
    
    label = paste0(
      source,
      " — ",
      threshold_display
    ),
    
    effect = OR_recalculated_exp_beta,
    lower  = OR_CI95_low_recalculated,
    upper  = OR_CI95_high_recalculated,
    
    estimate_text = sprintf(
      "%.2f (%.2f–%.2f)",
      effect,
      lower,
      upper
    ),
    
    p_text = sprintf(
      "%.3f",
      p_recalculated_from_beta_se
    )
  ) %>%
  arrange(
    taxon,
    source_order,
    threshold_order
  ) %>%
  select(
    source,
    taxon,
    analysis_role,
    instrument_threshold,
    nsnp_output,
    threshold_display,
    threshold_role,
    source_order,
    threshold_order,
    label,
    effect,
    lower,
    upper,
    estimate_text,
    p_text
  )

plotdata_path <- file.path(
  "06_figures",
  "08_plot_data",
  "Figure3_exposure_source_plot_data.csv"
)

readr::write_csv(
  fig3_dat,
  plotdata_path
)

stopifnot(file.exists(plotdata_path))

View(fig3_dat)
# ============================================================
# 3. Panel a: Veillonella
# ============================================================

fig3A_dat <- fig3_dat %>%
  filter(taxon == "Veillonella") %>%
  mutate(
    
    label_md = case_when(
      source == "MiBioGen" & threshold_order == 1 ~
        "MiBioGen — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "MiBioGen" & threshold_order == 2 ~
        "MiBioGen — *P* ≤ 1 × 10<sup>−5</sup>",
      
      source == "Swedish 2026" & threshold_order == 1 ~
        "Swedish 2026 — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "Swedish 2026" & threshold_order == 2 ~
        "Swedish 2026 — *P* ≤ 1 × 10<sup>−5</sup>",
      
      source == "DMP 2022" & threshold_order == 1 ~
        "DMP 2022 — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "DMP 2022" & threshold_order == 2 ~
        "DMP 2022 — *P* ≤ 1 × 10<sup>−5</sup>"
    ),
    
    # source groups之间留视觉间隔
    y_pos = case_when(
      source == "MiBioGen" & threshold_order == 1 ~ 6.0,
      source == "MiBioGen" & threshold_order == 2 ~ 5.0,
      
      source == "Swedish 2026" & threshold_order == 1 ~ 3.6,
      source == "Swedish 2026" & threshold_order == 2 ~ 2.6,
      
      source == "DMP 2022" & threshold_order == 1 ~ 1.2,
      source == "DMP 2022" & threshold_order == 2 ~ 0.2
    )
  )
#panel A 图
pA <- ggplot(
  fig3A_dat,
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
    breaks = c(
      0.6,
      0.75,
      1,
      1.25,
      1.5,
      1.75
    ),
    labels = c(
      "0.6",
      "0.75",
      "1",
      "1.25",
      "1.5",
      "1.75"
    )
  ) +
  coord_cartesian(
    xlim = c(0.55, 1.75)
  ) +
  scale_y_continuous(
    breaks = fig3A_dat$y_pos,
    labels = fig3A_dat$label_md,
    limits = c(-0.1, 6.5)
  ) +
  labs(
    title = "a. *Veillonella*",
    x = "Odds ratio for gastric cancer",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = ggtext::element_markdown(
      face = "bold",
      size = 13
    ),
    axis.text.y = ggtext::element_markdown(
      size = 9.5
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
# 4. Panel b: Veillonellaceae
# ============================================================

fig3B_dat <- fig3_dat %>%
  filter(taxon == "Veillonellaceae") %>%
  mutate(
    
    label_md = case_when(
      source == "MiBioGen" & threshold_order == 1 ~
        "MiBioGen — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "MiBioGen" & threshold_order == 2 ~
        "MiBioGen — *P* ≤ 1 × 10<sup>−5</sup>",
      
      source == "Swedish 2026" & threshold_order == 1 ~
        "Swedish 2026 — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "Swedish 2026" & threshold_order == 2 ~
        "Swedish 2026 — *P* ≤ 1 × 10<sup>−5</sup>",
      
      source == "DMP 2022" & threshold_order == 1 ~
        "DMP 2022 — *P* ≤ 5 × 10<sup>−6</sup>",
      
      source == "DMP 2022" & threshold_order == 2 ~
        "DMP 2022 — *P* ≤ 1 × 10<sup>−5</sup>"
    ),
    
    y_pos = case_when(
      source == "MiBioGen" & threshold_order == 1 ~ 6.0,
      source == "MiBioGen" & threshold_order == 2 ~ 5.0,
      
      source == "Swedish 2026" & threshold_order == 1 ~ 3.6,
      source == "Swedish 2026" & threshold_order == 2 ~ 2.6,
      
      source == "DMP 2022" & threshold_order == 1 ~ 1.2,
      source == "DMP 2022" & threshold_order == 2 ~ 0.2
    )
  )
#panel B 图
pB <- ggplot(
  fig3B_dat,
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
    breaks = c(
      0.6,
      0.75,
      1,
      1.25,
      1.5,
      1.75
    ),
    labels = c(
      "0.6",
      "0.75",
      "1",
      "1.25",
      "1.5",
      "1.75"
    )
  ) +
  coord_cartesian(
    xlim = c(0.55, 1.75)
  ) +
  scale_y_continuous(
    breaks = fig3B_dat$y_pos,
    labels = fig3B_dat$label_md,
    limits = c(-0.1, 6.5)
  ) +
  labs(
    title = "b. Veillonellaceae",
    x = "Odds ratio for gastric cancer",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.y = ggtext::element_markdown(
      size = 9.5
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
#先保存两张调试图
ggplot2::ggsave(
  filename = file.path(
    "06_figures",
    "09_outputs",
    "Figure3A_Veillonella.pdf"
  ),
  plot = pA,
  device = cairo_pdf,
  width = 7.2,
  height = 3.6,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    "06_figures",
    "09_outputs",
    "Figure3B_Veillonellaceae.pdf"
  ),
  plot = pB,
  device = cairo_pdf,
  width = 7.2,
  height = 3.6,
  units = "in"
)
#png
ggplot2::ggsave(
  filename = file.path(
    "06_figures",
    "09_outputs",
    "Figure3A_Veillonella.png"
  ),
  plot = pA,
  width = 7.2,
  height = 3.6,
  units = "in",
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    "06_figures",
    "09_outputs",
    "Figure3B_Veillonellaceae.png"
  ),
  plot = pB,
  width = 7.2,
  height = 3.6,
  units = "in",
  dpi = 300
)
#合图前调整：
# Panel a 去掉 x-axis title
pA <- pA +
  labs(x = NULL)

# Panel b 保留总 x-axis title
pB <- pB +
  labs(x = "Odds ratio for gastric cancer")
#合并
Figure3 <- pA / pB +
  patchwork::plot_layout(
    heights = c(1, 1)
  )

final_pdf <- file.path(
  "06_figures",
  "10_final",
  "Figure3_exposure_source_robustness.pdf"
)

ggplot2::ggsave(
  filename = final_pdf,
  plot = Figure3,
  device = cairo_pdf,
  width = 7.2,
  height = 6.4,
  units = "in"
)

final_png <- file.path(
  "06_figures",
  "10_final",
  "Figure3_exposure_source_robustness.png"
)

ggplot2::ggsave(
  filename = final_png,
  plot = Figure3,
  width = 7.2,
  height = 6.4,
  units = "in",
  dpi = 600
)
#检查
stopifnot(file.exists(final_pdf))
stopifnot(file.exists(final_png))

message(
  "Figure 3 script completed successfully. ",
  "Final PDF and PNG were generated."
)