# Reproduction instructions

This file explains how to reproduce the analysis from
*The State of Rental Discrimination* (Verhaeghe et al., forthcoming).

---

## Prerequisites

- **R** >= 4.3.0 (https://cran.r-project.org/)
- **RStudio** (optional but recommended; open `meta-analysis_housing.Rproj`)
- All R packages install automatically via `pak` blocks at the top of each script.
  An internet connection is required on first run.
- **Stage 5 only** (Bayesian models): install CmdStan after opening R:
  ```r
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  cmdstanr::install_cmdstan()
  ```
  Windows users also require Rtools: https://cran.r-project.org/bin/windows/Rtools/

---

## Data

The raw study-level workbook used in Stage 0 is not included in this repository
and is available from the corresponding author upon reasonable request.

The analysis-ready derived dataset (`1_data/public/housing_meta.rds` and
`1_data/public/housing_meta.csv`) is included and is sufficient to reproduce
**Stages 1-8** without the raw workbook.

---

## Stage order

Run scripts from the repository root, in order:

| Stage | Script | Input | Notes |
|---|---|---|---|
| 0 | `2_code/0_wrangling.R` | Raw workbook (not public) | Produces `housing_meta.rds`. Skip if starting from the derived dataset. |
| 1 | `2_code/1_descriptives.R` | `housing_meta.rds` | Table 1 and appendix country table |
| 2 | `2_code/2_meta_analysis.R` | `housing_meta.rds` | Pooled meta-analysis and PET/PEESE selection |
| 3 | `2_code/3_reg_remra.R` | `housing_meta.rds` | RE meta-regression |
| 4 | `2_code/4_reg_uwlsmra.R` | `housing_meta.rds` | UWLS meta-regression |
| 5 | `2_code/5_reg_bhmra.R` | `housing_meta.rds` | Bayesian hierarchical MRA -- **can take several hours** |
| 6 | `2_code/6_pub_bias.R` | `housing_meta.rds` | Publication bias diagnostics |
| 7 | `2_code/7_figures.R` | Stage outputs | Manuscript figures |
| 8 | `2_code/8_country_estimates.R` | `housing_meta.rds` | Country-level estimates (supplementary) |

From the command line:

```bash
Rscript 2_code/1_descriptives.R
Rscript 2_code/2_meta_analysis.R
Rscript 2_code/3_reg_remra.R
Rscript 2_code/4_reg_uwlsmra.R
Rscript 2_code/5_reg_bhmra.R          # optional; requires CmdStan
Rscript 2_code/6_pub_bias.R
Rscript 2_code/7_figures.R
Rscript 2_code/8_country_estimates.R  # optional supplementary
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
