Markdown
# Gastric Cancer–Gut Microbiome Mendelian Randomization

This repository contains the analysis code, derived results, publication-ready
plot data, figures, and tables for a Mendelian randomization (MR) study
investigating the potential causal relationship between gut microbiome traits
and gastric cancer.

The repository is organized to separate analysis scripts from derived
results and publication outputs.

---

## Study overview

This study investigates the potential causal relationship between selected
gut microbiome traits and gastric cancer using a Mendelian randomization
framework.

The core analyses include:

- Forward MR: gut microbiome traits → gastric cancer
- Reverse MR: gastric cancer → gut microbiome traits
- Threshold sensitivity analyses
- FinnGen outcome sensitivity analyses
- East Asian sensitivity analyses
- Exposure-source replication
- MR-Egger intercept testing
- Heterogeneity assessment
- MR-PRESSO analysis
- Leave-one-out analysis
- Single-SNP analysis

The principal microbiome exposures represented in the analysis include
Veillonella and Veillonellaceae.

---

## Repository structure

```text
gastric_cancer_gut_microbiome_MR/
│
├── README.md
├── .gitignore
│
├── 04_R_scripts/
│   ├── 00_functions.R
│   ├── 01_prepare_exposure.R
│   ├── 02_prepare_outcome.R
│   ├── 03_harmonise.R
│   ├── 04_primary_MR.R
│   ├── 05_threshold_sensitivity.R
│   ├── 06_finngen_sensitivity.R
│   ├── 07_east_asian_sensitivity.R
│   ├── 08_reverse_MR.R
│   ├── 09_summary_results.R
│   ├── 10_exploratory_exposure_replication.R
│   ├── 11_core_audit_checks.R
│   ├── 12_Table2_core_MR.R
│   ├── 13_Figure2_core_MR.R
│   ├── 14_Figure3_exposure_source_robustness.R
│   ├── 15_Table1_GWAS_characteristics&Analytical_roles.R
│   ├── 16_supplementary_tables_1to3.R
│   ├── 17_supplementary_tables_4to6.R
│   ├── 18_1_supplementary_figures2-4.R
│   └── 18_supplementary_figures.R
│
├── 05_results/
│   ├── 01_primary_MR/
│   ├── 02_threshold_sensitivity/
│   ├── 03_finngen_sensitivity/
│   ├── 04_east_asian_sensitivity/
│   ├── 05_reverse_MR/
│   ├── 06_summary/
│   └── 07_exposure_source_replication/
│
├── 06_figures/
│   ├── 08_plot_data/
│   └── 10_final/
│
└── 07_tables/
    ├── 02_publication_data/
    └── 03_outputs/

Analysis workflow
The analysis follows the general workflow:
GWAS exposure data
        │
        ▼
01_prepare_exposure.R
        │
        ▼
02_prepare_outcome.R
        │
        ▼
03_harmonise.R
        │
        ▼
04_primary_MR.R
        │
        ├── 05_threshold_sensitivity.R
        ├── 06_finngen_sensitivity.R
        ├── 07_east_asian_sensitivity.R
        ├── 08_reverse_MR.R
        └── 10_exploratory_exposure_replication.R
                │
                ▼
        11_core_audit_checks.R
                │
                ▼
        Summary results
                │
        ┌───────┴────────┐
        ▼                ▼
    Figures            Tables
Additional scripts generate the main manuscript tables and figures and
supplementary outputs.

Analysis scripts

Core data preparation
| Script | Purpose |
|---|---|
| `00_functions.R` | Shared functions used throughout the analysis |
| `01_prepare_exposure.R` | Preparation of microbiome exposure instruments |
| `02_prepare_outcome.R` | Preparation of gastric cancer outcome data |
| `03_harmonise.R` | Harmonisation of exposure and outcome datasets |

Main MR analysis
| Script | Purpose |
|---|---|
| `04_primary_MR.R` | Primary MR analysis |
| `05_threshold_sensitivity.R` | Alternative instrument-selection threshold analysis |
| `06_finngen_sensitivity.R` | FinnGen outcome sensitivity analysis |
| `07_east_asian_sensitivity.R` | East Asian sensitivity analysis |
| `08_reverse_MR.R` | Reverse-direction MR analysis |
| `09_summary_results.R` | Summary of MR results |
| `10_exploratory_exposure_replication.R` | Exposure-source replication analysis |

Audit and publication outputs
| Script | Purpose |
|---|---|
| `11_core_audit_checks.R` | Core analytical and result consistency checks |
| `12_Table2_core_MR.R` | Generation of the core MR table |
| `13_Figure2_core_MR.R` | Generation of the core bidirectional MR figure |
| `14_Figure3_exposure_source_robustness.R` | Exposure-source robustness analysis figure |
| `15_Table1_GWAS_characteristics&Analytical_roles.R` | GWAS characteristics and analytical roles |
| `16_supplementary_tables_1to3.R` | Supplementary Tables 1–3 |
| `17_supplementary_tables_4to6.R` | Supplementary Tables 4–6 |
| `18_1_supplementary_figures2-4.R` | Supplementary Figures 2–4 |
| `18_supplementary_figures.R` | Supplementary figure generation |

Results
05_results/ contains derived analytical results organized according to
the analysis stage.
Primary MR
05_results/01_primary_MR/
Contains the primary forward MR results and corresponding diagnostic outputs,
including:
- MR estimates
- MR-Egger intercepts
- heterogeneity statistics
- MR-PRESSO results
- leave-one-out results
- single-SNP results
Threshold sensitivity
05_results/02_threshold_sensitivity/
Contains results obtained using the alternative instrument-selection
threshold.
FinnGen sensitivity
05_results/03_finngen_sensitivity/
Contains sensitivity analyses using FinnGen gastric cancer outcome data.
East Asian sensitivity
05_results/04_east_asian_sensitivity/
Contains analyses using the available East Asian outcome dataset.
Reverse MR
05_results/05_reverse_MR/
Contains reverse-direction MR analyses evaluating: gastric cancer → gut microbiome

Summary
05_results/06_summary/
Contains consolidated MR and QC summary files.
Exposure-source replication
05_results/07_exposure_source_replication/
Contains analyses evaluating the robustness of the results to an independent
or alternative exposure GWAS source.

Figures
06_figures/ contains plot data and final figures.
Plot data
06_figures/08_plot_data/
Contains the data used to generate selected manuscript figures.
Final figures
06_figures/10_final/
Contains final manuscript figures and supplementary figures in publication
formats such as PDF and PNG.

Tables
07_tables/ contains publication-oriented tables.
Publication data
07_tables/02_publication_data/
Contains machine-readable CSV files corresponding to the main and
supplementary publication tables.
Table outputs
07_tables/03_outputs/
Contains formatted Word versions of the main and supplementary tables.

Data availability
Raw GWAS datasets and other large or externally sourced input datasets are
not included in this repository.
The repository instead contains the analysis scripts and selected derived
outputs required to document the analytical workflow and reproduce the
reported downstream results from the corresponding source data.
Users should obtain the original GWAS datasets from their respective public
repositories or study-specific data-access mechanisms.

Files intentionally excluded from the public repository
The following materials are intentionally excluded:
- raw GWAS data
- private or intermediate datasets
- administrative and registration materials
- manuscript working files
- historical project snapshots
- private reference libraries
- local R/RStudio files
- temporary files
These exclusions are enforced through .gitignore.

Reproducibility
The analysis scripts are retained in 04_R_scripts/.
Derived results, figure data, figures, and publication tables are provided in
their corresponding output directories.
The repository therefore separates:
Analysis code
      ↓
Derived results
      ↓
Plot data
      ↓
Final figures
      ↓
Publication tables
This structure is intended to make the relationship between analytical
scripts and reported outputs transparent.

Software
The analyses were conducted using R and relevant R packages for Mendelian
randomization, GWAS processing, statistical analysis, and visualization.
Exact package versions and computational-environment records are maintained
in the project materials outside the public repository where appropriate.

Citation
If you use this repository or its analysis workflow, please cite the
corresponding study when the manuscript is published.
