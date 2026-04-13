# Housing Discrimination Meta-Analysis

Replication code and analysis-ready data for:

> Verhaeghe, P. P., Devos, L., Ghekiere, A., Baert, S., & Lippens, L. (*forthcoming*).
> The state of rental discrimination: A meta-analytical comparison of five discrimination
> grounds in the housing market.
> OSF: [https://osf.io/t7hbe](https://osf.io/t7hbe)

---

## Repository contents

| Path | Description |
|---|---|
| `1_data/public/housing_meta.csv` | Analysis-ready dataset (475 effect sizes, 114 studies) |
| `1_data/public/housing_meta.rds` | Same dataset in R binary format (loaded by all downstream scripts) |
| `2_code/1_wrangling.R` | Stage 1 — data intake, effect-size computation, sample flags |
| `2_code/2_descriptives.R` | Stage 2 — Table 1 and appendix country table |
| `2_code/3_meta_analysis.R` | Stage 3 — pooled meta-analysis and PET/PEESE selection |
| `2_code/4_reg_remra.R` | Stage 4 — random-effects meta-regression (RE-MRA) |
| `2_code/5_reg_uwlsmra.R` | Stage 5 — UWLS meta-regression with bias correction |
| `2_code/6_reg_bhmra.R` | Stage 6 — Bayesian hierarchical meta-regression (BH-MRA) |
| `2_code/7_pub_bias.R` | Stage 7 — publication bias diagnostics |
| `2_code/8_figures.R` | Stage 8 — manuscript figure export |
| `2_code/9_country_estimates.R` | Stage 9 — country-level estimates (supplementary) |
| `2_code/functions/` | Shared helper functions (core, reg, io) |
| `REPRODUCE.md` | Step-by-step reproduction instructions |
| `CLAUDE.md` | Detailed project reference (data dictionary, pipeline, conventions) |
| `AGENTS.md` | Analysis agent instructions and coding standards |

The raw study-level workbook (`1_data/260206_housing.xlsx`) is not included in this
repository. It is available from the corresponding author upon reasonable request.

---

## How to reproduce

See [`REPRODUCE.md`](REPRODUCE.md) for step-by-step instructions, including
package installation, `cmdstanr` setup (required for Stage 6), and the recommended
run order.

**Quick start** (Stages 2-9 only, skipping Stage 1 which requires the raw workbook):

```r
# All packages install automatically via pak blocks at the top of each script.
# Run scripts in order from your R console or RStudio:
source("2_code/2_descriptives.R")
source("2_code/3_meta_analysis.R")
# ... etc.
```

---

## Data

`housing_meta.csv` / `housing_meta.rds` is the canonical analysis-ready dataset
derived from `1_data/260206_housing.xlsx` by `2_code/1_wrangling.R`. It contains
475 effect sizes (positive response log-ratios) from 114 correspondence audit studies
across five discrimination grounds (ethno-racial origin, gender, health and disability,
sexual orientation, social origin).

The full replication archive — including final figures and tables — is available on OSF:
[https://osf.io/t7hbe](https://osf.io/t7hbe)

---

## Software

- **R** >= 4.3.0
- Key packages: `metafor`, `fixest`, `brms`, `cmdstanr`, `modelsummary`, `ggplot2`
- `cmdstanr` requires a separate CmdStan installation: `cmdstanr::install_cmdstan()`
- All other packages install automatically via `pak` blocks at the top of each script.

---

## License

Code: [MIT License](LICENSE)
Data (`housing_meta.csv`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
