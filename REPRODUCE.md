# Reproduction instructions

This file explains how to reproduce the analysis from
*The State of Rental Discrimination* (Verhaeghe et al., forthcoming).

---

## Prerequisites

- **R** >= 4.3.0 (https://cran.r-project.org/)
- **RStudio** (optional but recommended; open `meta-analysis_housing.Rproj`)
- All R packages install automatically via `pak` blocks at the top of each script.
  An internet connection is required on first run.
- **Stage 6 only** (Bayesian models): install CmdStan after opening R:
  ```r
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  cmdstanr::install_cmdstan()
  ```
  Windows users also require Rtools: https://cran.r-project.org/bin/windows/Rtools/

---

## Data

The raw study-level workbook (`1_data/260206_housing.xlsx`) is not included in this
repository and is available from the corresponding author upon reasonable request.

The analysis-ready derived dataset (`1_data/public/housing_meta.rds` and
`1_data/public/housing_meta.csv`) is included and is sufficient to reproduce
**Stages 2-9** without the raw workbook.

---

## Stage order

Run scripts from the repository root, in order:

| Stage | Script | Input | Notes |
|---|---|---|---|
| 1 | `2_code/1_wrangling.R` | `1_data/260206_housing.xlsx` | Requires raw data (not public). Produces `housing_meta.rds`. Skip if starting from the derived dataset. |
| 2 | `2_code/2_descriptives.R` | `housing_meta.rds` | Table 1 and appendix country table |
| 3 | `2_code/3_meta_analysis.R` | `housing_meta.rds` | Pooled meta-analysis and PET/PEESE selection |
| 4 | `2_code/4_reg_remra.R` | `housing_meta.rds` | RE meta-regression |
| 5 | `2_code/5_reg_uwlsmra.R` | `housing_meta.rds` | UWLS meta-regression |
| 6 | `2_code/6_reg_bhmra.R` | `housing_meta.rds` | Bayesian hierarchical MRA -- **can take several hours** |
| 7 | `2_code/7_pub_bias.R` | `housing_meta.rds` | Publication bias diagnostics |
| 8 | `2_code/8_figures.R` | Stage outputs | Manuscript figures |
| 9 | `2_code/9_country_estimates.R` | Stage 4/5 outputs | Country-level estimates (supplementary) |

From the command line:

```bash
Rscript 2_code/2_descriptives.R
Rscript 2_code/3_meta_analysis.R
Rscript 2_code/4_reg_remra.R
Rscript 2_code/5_reg_uwlsmra.R
Rscript 2_code/6_reg_bhmra.R          # optional; requires CmdStan
Rscript 2_code/7_pub_bias.R
Rscript 2_code/8_figures.R
Rscript 2_code/9_country_estimates.R  # optional supplementary
```

---

## Expected outputs

All outputs are written to `3_figures/` (figures) and `4_tables/` (tables and model objects).
A frozen snapshot of the final outputs is archived on OSF: https://osf.io/t7hbe

---

## Seed

All stages use `set.seed(8888)`.

---

## CSV export

Outputs are saved as `.rds` by default. To also export `.csv` files, set before running:

```r
options(housing.export.csv = TRUE)
```

or set environment variable `HOUSING_EXPORT_CSV=true`.
