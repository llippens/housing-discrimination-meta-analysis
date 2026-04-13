# The state of rental discrimination

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
| `2_code/1_descriptives.R` | Stage 1 — Table 1 and appendix country table |
| `2_code/2_meta_analysis.R` | Stage 2 — pooled meta-analysis and PET/PEESE selection |
| `2_code/3_reg_remra.R` | Stage 3 — random-effects meta-regression (RE-MRA) |
| `2_code/4_reg_uwlsmra.R` | Stage 4 — UWLS meta-regression with bias correction |
| `2_code/5_reg_bhmra.R` | Stage 5 — Bayesian hierarchical meta-regression (BH-MRA) |
| `2_code/6_pub_bias.R` | Stage 6 — publication bias diagnostics |
| `2_code/7_figures.R` | Stage 7 — manuscript figure export |
| `2_code/8_country_estimates.R` | Stage 8 — country-level estimates (supplementary) |
| `2_code/functions/` | Shared helper functions (core, reg, io) |
| `REPRODUCE.md` | Step-by-step reproduction instructions |

The public data consist of cleaned and wrangled metadata; the original workbook is not included in this
repository. It is available from the corresponding author upon request.

---

## How to reproduce

See [`REPRODUCE.md`](REPRODUCE.md) for step-by-step instructions, including
package installation, `cmdstanr` setup (required for Stage 5), and the recommended
run order.

**Quick start** (Stages 1-8, starting from the derived dataset):

```r
# All packages install automatically via pak blocks at the top of each script.
# Run scripts in order from your R console or RStudio:
source("2_code/1_descriptives.R")
source("2_code/2_meta_analysis.R")
# ... etc.
```

---

## Data

`housing_meta.csv` / `housing_meta.rds` is the canonical analysis-ready dataset
containing 475 effect sizes (positive response log-ratios) from 114 correspondence
audit studies across five discrimination grounds (ethno-racial origin, gender,
health and disability, sexual orientation, social origin).

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
