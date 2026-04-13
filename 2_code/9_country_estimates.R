# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "dplyr", "purrr", "tibble",
  "metafor", "clubSandwich", "fixest", "marginaleffects",
  "writexl"
)

missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_models.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_design_utils.R"))
source(file.path(here::here(), "2_code", "functions", "reg", "f_metareg_helpers.R"))
source(file.path(here::here(), "2_code", "functions", "reg", "f_uwls_helpers.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root    <- here::here()
dir.public  <- file.path(dir.root, "1_data", "public")
dir.tables  <- file.path(dir.root, "4_tables")
dir.cty.out <- file.path(dir.tables, "9_country_estimates")
ensure_dir(dir.cty.out)

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/1_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "prratio_ln", "prratio_ln_v", "study", "cluster", "ground_abbr",
  "country", "treatment_group", "year_mid", "info_binary", "agent_type",
  "matched_design", "audit_type", "callback_type", "airbnb",
  "can_gender", "can_education", "can_employment",
  "can_migrant_generation", "peer_reviewed", "language"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Model setup ####
housing.meta.by.ground    <- split(housing.meta, as.character(housing.meta$ground_abbr))
housing.meta.scope.data   <- housing.meta.by.ground

factor.cols <- c(
  "country", "info_binary", "agent_type", "matched_design", "audit_type",
  "treatment_group", "callback_type", "airbnb", "can_gender",
  "can_education", "can_employment", "can_migrant_generation",
  "peer_reviewed", "language"
)
numeric.cols <- c("year_mid")

# Helper functions ####

cty_reg_candidates_for_scope <- function(scope_name) {
  terms <- c(
    "year_mid", "info_binary", "agent_type",
    "matched_design", "audit_type", "callback_type",
    "airbnb", "can_education", "can_employment",
    "peer_reviewed", "language"
  )
  if (scope_name != "gen") terms <- c(terms, "can_gender")
  if (scope_name == "ero")  terms <- c(terms, "can_migrant_generation")
  unique(terms)
}

# Observed estimates ####
# Inverse-variance weighted mean per country (plain pooling, no model)
compute_obs_estimates <- function(data, scope_name) {
  data %>%
    dplyr::filter(
      is.finite(prratio_ln),
      is.finite(prratio_ln_v),
      prratio_ln_v > 0,
      !is.na(country),
      as.character(country) != ""
    ) %>%
    dplyr::group_by(country = as.character(country)) %>%
    dplyr::summarise(
      ground_abbr    = scope_name,
      k              = dplyr::n(),
      n_studies      = dplyr::n_distinct(study),
      obs_prratio_ln = sum(prratio_ln / prratio_ln_v) / sum(1 / prratio_ln_v),
      obs_se         = sqrt(1 / sum(1 / prratio_ln_v)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      obs_ci_low  = obs_prratio_ln - 1.96 * obs_se,
      obs_ci_high = obs_prratio_ln + 1.96 * obs_se
    )
}

# RE-MRA helpers ####

# predict.rma.mv: standard vcov, not CR2 — acceptable for supplementary predictions.
build_remra_newmods <- function(fit, data_use, country_val) {
  x_cols   <- colnames(fit$X)
  mod_cols <- x_cols[x_cols != "intrcpt"]

  newmods_v <- setNames(rep(0, length(mod_cols)), mod_cols)

  cty_col <- paste0("country", country_val)
  if (cty_col %in% mod_cols) newmods_v[cty_col] <- 1L

  for (col in mod_cols) {
    if (col %in% names(data_use) && is.numeric(data_use[[col]])) {
      valid <- is.finite(data_use[[col]]) & is.finite(data_use$prratio_ln_v) & data_use$prratio_ln_v > 0
      if (any(valid)) {
        w <- 1 / data_use$prratio_ln_v[valid]
        newmods_v[col] <- sum(data_use[[col]][valid] * w) / sum(w)
      }
    }
  }

  matrix(newmods_v, nrow = 1L, dimnames = list(NULL, mod_cols))
}

# Fit RE-MRA with country as geography moderator; return per-country predictions.
fit_remra_country <- function(data, scope_name, factor_cols, numeric_cols, with_covs) {
  col_prefix <- if (with_covs) "remra_cov" else "remra_cty"

  data_use <- prepare_metareg_data(data, factor_cols, numeric_cols)
  # split() preserves all parent factor levels; drop empty ones to avoid collinearity warnings.
  data_use <- droplevels(data_use)
  if ("country" %in% names(data_use) && is.factor(data_use$country) &&
      "United States" %in% levels(data_use$country)) {
    data_use$country <- relevel(data_use$country, ref = "United States")
  }

  n_studies    <- dplyr::n_distinct(data_use$study[!is.na(data_use$study)])
  random_terms <- infer_random_terms(data_use)
  countries    <- sort(unique(as.character(data_use$country[!is.na(data_use$country)])))

  make_na_tbl <- function(status_val) {
    tibble::tibble(country = countries) %>%
      dplyr::mutate(
        !!paste0(col_prefix, "_est")     := NA_real_,
        !!paste0(col_prefix, "_se")      := NA_real_,
        !!paste0(col_prefix, "_ci_low")  := NA_real_,
        !!paste0(col_prefix, "_ci_high") := NA_real_,
        !!paste0(col_prefix, "_conv_ok") := NA,
        !!paste0(col_prefix, "_status")  := status_val
      )
  }

  fixed_terms <- if (with_covs) {
    cov_candidates <- cty_reg_candidates_for_scope(scope_name)
    cov_candidates <- cov_candidates[cov_candidates %in% names(data_use)]
    cov_use <- cov_candidates[vapply(cov_candidates, function(v) is_informative_fixed(data_use[[v]]), logical(1))]
    c("country", cov_use)
  } else {
    "country"
  }

  df_stats <- calc_design_df(data_use, fixed_terms, random_terms, random_df_method = "study_only")
  if (df_stats$total_df > n_studies) return(make_na_tbl("df_exceeded"))

  formula_mod <- stats::as.formula(paste("~", paste(fixed_terms, collapse = " + ")))
  fit <- tryCatch(fit_rma_mv_mod(data_use, mods = formula_mod), error = function(e) e)
  if (inherits(fit, "error")) return(make_na_tbl("error"))

  # fit$code: nlminb returns 0 for success; treat NULL (some metafor builds / optimizers) as converged
  code_ok <- is.null(fit$code) || isTRUE(fit$code == 0)
  conv_ok <- code_ok && !any(is.na(fit$beta)) && all(is.finite(fit$beta))

  purrr::map_dfr(countries, function(ctr) {
    newmods_row <- tryCatch(build_remra_newmods(fit, data_use, ctr), error = function(e) NULL)
    pred <- if (!is.null(newmods_row)) {
      tryCatch(predict(fit, newmods = newmods_row), error = function(e) NULL)
    } else {
      NULL
    }

    row_status <- dplyr::case_when(
      !conv_ok         ~ "no_conv",
      is.null(pred)    ~ "prediction_error",
      TRUE             ~ "ok"
    )

    tibble::tibble(country = ctr) %>%
      dplyr::mutate(
        !!paste0(col_prefix, "_est")     := if (!is.null(pred)) as.numeric(pred$pred)  else NA_real_,
        !!paste0(col_prefix, "_se")      := if (!is.null(pred)) as.numeric(pred$se)    else NA_real_,
        !!paste0(col_prefix, "_ci_low")  := if (!is.null(pred)) as.numeric(pred$ci.lb) else NA_real_,
        !!paste0(col_prefix, "_ci_high") := if (!is.null(pred)) as.numeric(pred$ci.ub) else NA_real_,
        !!paste0(col_prefix, "_conv_ok") := conv_ok,
        !!paste0(col_prefix, "_status")  := row_status
      )
  })
}

# UWLS helpers ####

# Fit UWLS with country in RHS (not absorbed as FE); treatment_group stays as FE.
fit_uwls_country <- function(data, scope_name, model_type, factor_cols, numeric_cols, with_covs) {
  pet_term <- if (model_type == "pet") "prratio_ln_se" else "prratio_ln_v"
  data_use <- prepare_uwls_data(data, factor_cols, numeric_cols)
  data_use <- droplevels(data_use)

  if ("country" %in% names(data_use) && is.factor(data_use$country) &&
      "United States" %in% levels(data_use$country)) {
    data_use$country <- relevel(data_use$country, ref = "United States")
  }

  n_rows    <- nrow(data_use)
  n_studies <- dplyr::n_distinct(data_use$study[!is.na(data_use$study)])
  fe_terms  <- "treatment_group"
  fe_use    <- fe_terms[
    fe_terms %in% names(data_use) &
      vapply(fe_terms, function(v) is_informative_fixed(data_use[[v]]), logical(1))
  ]

  reg_use <- if (with_covs) {
    cov_candidates <- cty_reg_candidates_for_scope(scope_name)
    cov_candidates <- cov_candidates[cov_candidates %in% names(data_use)]
    cov_use <- cov_candidates[vapply(cov_candidates, function(v) is_informative_fixed(data_use[[v]]), logical(1))]
    setdiff(c("country", cov_use), fe_use)
  } else {
    setdiff("country", fe_use)
  }

  if (n_rows < 2) {
    return(list(model = NULL, data = data_use, status = "insufficient_n", conv_ok = FALSE, pet_p = NA_real_))
  }
  if (n_studies <= 1) {
    return(list(model = NULL, data = data_use, status = "insufficient_random_levels", conv_ok = FALSE, pet_p = NA_real_))
  }

  df_stats <- calc_uwls_df(data_use, reg_terms = reg_use, fe_terms = fe_use, pet_term = pet_term)
  if (df_stats$total_df > n_studies) {
    return(list(model = NULL, data = data_use, status = "df_exceeded", conv_ok = FALSE, pet_p = NA_real_))
  }

  formula_mod <- build_uwls_formula(
    reg_terms         = reg_use,
    interaction_terms = character(0),
    factor_cols       = factor_cols,
    pet_term          = pet_term,
    fe_terms          = fe_use
  )

  fit <- tryCatch(
    fixest::feols(
      fml     = formula_mod,
      data    = data_use,
      weights = data_use$uwls_weight,
      vcov    = fixest::vcov_cluster("study")
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(model = NULL, data = data_use, status = "error", conv_ok = FALSE, pet_p = NA_real_))
  }

  pet_p <- NA_real_
  if (model_type == "pet") {
    ct <- tryCatch(fixest::coeftable(fit, vcov = ~study, ssc = ssc_strict), error = function(e) NULL)
    if (!is.null(ct)) {
      pet_row <- rownames(ct) == pet_term
      if (any(pet_row)) pet_p <- unname(ct[pet_row, 4])
    }
  }

  list(model = fit, data = data_use, status = "ok", conv_ok = TRUE, pet_p = pet_p)
}

# Predictions per country. Precision term fixed at prediction time to remove
# the publication-bias component (PET: se = 0.1; PEESE: v = 0.01).
# use_counterfactual = TRUE  (default): counterfactual marginal standardisation —
#   datagrid(..., grid_type = "counterfactual") creates n_obs rows with country
#   fixed, then avg_predictions averages weighted by 1/precision. Used for cty.
# use_counterfactual = FALSE: prediction at covariate means/modes — datagrid()
#   produces one representative row (mean of numerics, mode of factors). Avoids
#   delta-method SE inflation from cluster-robust vcov propagated over n_obs
#   correlated rows, which renders CIs uninformative with many covariates (cov).
get_uwls_country_preds <- function(fit_result, countries, model_type,
                                   use_counterfactual = TRUE) {
  if (is.null(fit_result$model)) {
    return(tibble::tibble(
      country = countries,
      est = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_,
      conv_ok = NA, status = fit_result$status
    ))
  }

  # Collinear country dummies predict at US reference level — silent and misleading; return NA.
  collin_vars <- fit_result$model$collin.var
  collin_ctrs <- if (length(collin_vars) > 0) {
    sub("^country", "", collin_vars[startsWith(collin_vars, "country")])
  } else {
    character(0)
  }

  purrr::map_dfr(countries, function(ctr) {
    if (ctr %in% collin_ctrs) {
      return(tibble::tibble(
        country = ctr,
        est = NA_real_, se = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_,
        conv_ok = FALSE, status = "collinear_dropped"
      ))
    }

    tryCatch({
      nd_args <- list(model = fit_result$model, country = ctr)
      if (model_type == "pet")   nd_args$prratio_ln_se <- 1e-1
      if (model_type == "peese") nd_args$prratio_ln_v  <- 1e-2
      if (use_counterfactual)    nd_args$grid_type      <- "counterfactual"

      nd <- do.call(marginaleffects::datagrid, nd_args)

      if (use_counterfactual) {
        # feols only stores formula columns; datagrid adds exactly one precision column.
        wts_col <- if (model_type == "pet") "prratio_ln_se" else "prratio_ln_v"
        pred <- marginaleffects::avg_predictions(fit_result$model, newdata = nd,
                                                 wts = 1 / nd[[wts_col]])
      } else {
        # Single representative row — no weighting needed.
        pred <- marginaleffects::avg_predictions(fit_result$model, newdata = nd)
      }
      pred <- pred[1L, ]

      tibble::tibble(
        country = ctr,
        est     = pred$estimate,
        se      = pred$std.error,
        ci_low  = pred$conf.low,
        ci_high = pred$conf.high,
        conv_ok = is.finite(pred$estimate) && is.finite(pred$std.error),
        status  = "ok"
      )
    }, error = function(e) {
      # datagrid errors on country levels absent from training data.
      tibble::tibble(
        country = ctr,
        est = NA_real_, se = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_,
        conv_ok = FALSE, status = "prediction_error"
      )
    })
  })
}

# PET/PEESE selection mirrors stage-5 rule (p < 0.10 → PEESE, else PET).
select_pet_peese <- function(pet_obj, peese_obj) {
  if (pet_obj$status == "ok" && is.finite(pet_obj$pet_p) && pet_obj$pet_p < 0.10) {
    list(selected = "peese", obj = peese_obj)
  } else if (pet_obj$status == "ok") {
    list(selected = "pet", obj = pet_obj)
  } else if (peese_obj$status == "ok") {
    list(selected = "peese", obj = peese_obj)
  } else {
    list(
      selected = "none",
      obj = list(model = NULL, data = NULL, status = pet_obj$status, conv_ok = FALSE, pet_p = NA_real_)
    )
  }
}

# Main execution ####
housing.meta.cty.list <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    message("Processing scope: ", scope_name)

    obs_tbl   <- compute_obs_estimates(df, scope_name)
    countries <- sort(obs_tbl$country)

    remra_cty_tbl <- fit_remra_country(df, scope_name, factor.cols, numeric.cols, with_covs = FALSE)
    remra_cov_tbl <- fit_remra_country(df, scope_name, factor.cols, numeric.cols, with_covs = TRUE)

    uwls_cty_pet   <- fit_uwls_country(df, scope_name, "pet",   factor.cols, numeric.cols, with_covs = FALSE)
    uwls_cty_peese <- fit_uwls_country(df, scope_name, "peese", factor.cols, numeric.cols, with_covs = FALSE)
    uwls_cov_pet   <- fit_uwls_country(df, scope_name, "pet",   factor.cols, numeric.cols, with_covs = TRUE)
    uwls_cov_peese <- fit_uwls_country(df, scope_name, "peese", factor.cols, numeric.cols, with_covs = TRUE)

    cty_sel <- select_pet_peese(uwls_cty_pet, uwls_cty_peese)
    cov_sel <- select_pet_peese(uwls_cov_pet, uwls_cov_peese)

    uwls_cty_preds <- get_uwls_country_preds(cty_sel$obj, countries, cty_sel$selected)
    uwls_cov_preds <- get_uwls_country_preds(cov_sel$obj, countries, cov_sel$selected,
                                             use_counterfactual = FALSE)

    uwls_cty_tbl <- uwls_cty_preds %>%
      dplyr::rename(
        uwls_cty_est     = est,     uwls_cty_se       = se,
        uwls_cty_ci_low  = ci_low,  uwls_cty_ci_high  = ci_high,
        uwls_cty_conv_ok = conv_ok, uwls_cty_status   = status
      ) %>%
      dplyr::mutate(uwls_cty_model_type = cty_sel$selected)

    uwls_cov_tbl <- uwls_cov_preds %>%
      dplyr::rename(
        uwls_cov_est     = est,     uwls_cov_se       = se,
        uwls_cov_ci_low  = ci_low,  uwls_cov_ci_high  = ci_high,
        uwls_cov_conv_ok = conv_ok, uwls_cov_status   = status
      ) %>%
      dplyr::mutate(uwls_cov_model_type = cov_sel$selected)

    obs_tbl %>%
      dplyr::left_join(remra_cty_tbl, by = "country") %>%
      dplyr::left_join(remra_cov_tbl, by = "country") %>%
      dplyr::left_join(uwls_cty_tbl,  by = "country") %>%
      dplyr::left_join(uwls_cov_tbl,  by = "country")
  }
)

# Combine scopes and back-transform ####
housing.meta.cty.estimates <- dplyr::bind_rows(housing.meta.cty.list) %>%
  dplyr::mutate(
    obs_prratio          = exp(obs_prratio_ln),
    obs_prratio_ci_low   = exp(obs_ci_low),
    obs_prratio_ci_high  = exp(obs_ci_high),
    remra_cty_prratio         = exp(remra_cty_est),
    remra_cty_prratio_ci_low  = exp(remra_cty_ci_low),
    remra_cty_prratio_ci_high = exp(remra_cty_ci_high),
    remra_cov_prratio         = exp(remra_cov_est),
    remra_cov_prratio_ci_low  = exp(remra_cov_ci_low),
    remra_cov_prratio_ci_high = exp(remra_cov_ci_high),
    uwls_cty_prratio         = exp(uwls_cty_est),
    uwls_cty_prratio_ci_low  = exp(uwls_cty_ci_low),
    uwls_cty_prratio_ci_high = exp(uwls_cty_ci_high),
    uwls_cov_prratio         = exp(uwls_cov_est),
    uwls_cov_prratio_ci_low  = exp(uwls_cov_ci_low),
    uwls_cov_prratio_ci_high = exp(uwls_cov_ci_high)
  ) %>%
  dplyr::select(
    ground_abbr, country, k, n_studies,
    obs_prratio_ln, obs_se, obs_ci_low, obs_ci_high,
    obs_prratio, obs_prratio_ci_low, obs_prratio_ci_high,
    remra_cty_est, remra_cty_se, remra_cty_ci_low, remra_cty_ci_high,
    remra_cty_conv_ok, remra_cty_status,
    remra_cty_prratio, remra_cty_prratio_ci_low, remra_cty_prratio_ci_high,
    remra_cov_est, remra_cov_se, remra_cov_ci_low, remra_cov_ci_high,
    remra_cov_conv_ok, remra_cov_status,
    remra_cov_prratio, remra_cov_prratio_ci_low, remra_cov_prratio_ci_high,
    uwls_cty_est, uwls_cty_se, uwls_cty_ci_low, uwls_cty_ci_high,
    uwls_cty_conv_ok, uwls_cty_status, uwls_cty_model_type,
    uwls_cty_prratio, uwls_cty_prratio_ci_low, uwls_cty_prratio_ci_high,
    uwls_cov_est, uwls_cov_se, uwls_cov_ci_low, uwls_cov_ci_high,
    uwls_cov_conv_ok, uwls_cov_status, uwls_cov_model_type,
    uwls_cov_prratio, uwls_cov_prratio_ci_low, uwls_cov_prratio_ci_high
  )

# Data dictionary ####

housing.meta.cty.dict <- tibble::tribble(
  ~column,                    ~section,             ~scale,        ~description,                                                                                                   ~notes,

  "ground_abbr",              "Identifiers",        "categorical", "Discrimination ground abbreviation.",                                                                          "ero = ethno-racial origin; gen = gender; hed = health and disability; seo = sexual orientation; soc = social origin.",
  "country",                  "Identifiers",        "categorical", "Country where the correspondence audit studies were conducted.",                                               "Country names follow the labels in housing_meta.rds. Reference country for RE-MRA and UWLS models is United States.",
  "k",                        "Identifiers",        "integer",     "Number of individual effect sizes (audit observations) in this country x ground cell.",                        "Each effect size corresponds to one audit pair or matched set. One study can contribute multiple effect sizes.",
  "n_studies",                "Identifiers",        "integer",     "Number of distinct audit studies contributing to this cell.",                                                  "Used as the effective sample size for the degrees-of-freedom gate that decides whether covariate models are attempted.",

  "obs_prratio_ln",           "Observed (IVW)",     "log PR ratio","Inverse-variance weighted (IVW) mean of the log positive-response-rate ratio across all effect sizes in the cell.", "Formula: sum(y_i / v_i) / sum(1 / v_i), where y_i is the log PR ratio and v_i is its sampling variance. No modelling, no bias correction, no clustering correction.",
  "obs_se",                   "Observed (IVW)",     "log PR ratio","Standard error of the IVW mean.",                                                                             "Formula: sqrt(1 / sum(1 / v_i)). Treats all effect sizes as independent; does not account for within-study clustering.",
  "obs_ci_low",               "Observed (IVW)",     "log PR ratio","Lower bound of the 95% confidence interval for the IVW mean (obs_prratio_ln - 1.96 * obs_se).",              NA_character_,
  "obs_ci_high",              "Observed (IVW)",     "log PR ratio","Upper bound of the 95% confidence interval for the IVW mean (obs_prratio_ln + 1.96 * obs_se).",              NA_character_,
  "obs_prratio",              "Observed (IVW)",     "PR ratio",    "IVW mean PR ratio on the original ratio scale: exp(obs_prratio_ln).",                                         "Values below 1 indicate net discrimination against the minority group in this country.",
  "obs_prratio_ci_low",       "Observed (IVW)",     "PR ratio",    "Lower bound of the 95% CI for the PR ratio: exp(obs_ci_low).",                                               NA_character_,
  "obs_prratio_ci_high",      "Observed (IVW)",     "PR ratio",    "Upper bound of the 95% CI for the PR ratio: exp(obs_ci_high).",                                              NA_character_,

  "remra_cty_est",            "RE-MRA: country only","log PR ratio","Random-effects meta-regression (RE-MRA) estimate for the country (log PR ratio). Model: rma.mv with ~ country as the only fixed moderator and nested random effects (study / cluster). Reference country = United States.", "Accounts for study-level clustering via random effects. Does NOT correct for publication bias. Prediction is the model-estimated log PR ratio for this country, all else equal.",
  "remra_cty_se",             "RE-MRA: country only","log PR ratio","Standard error of the RE-MRA country estimate (log scale).",                                               "Derived from predict.rma.mv using the standard (not CR2-adjusted) variance-covariance matrix. Appropriate for this supplementary prediction context.",
  "remra_cty_ci_low",         "RE-MRA: country only","log PR ratio","Lower bound of the 95% CI from predict.rma.mv (log scale).",                                               NA_character_,
  "remra_cty_ci_high",        "RE-MRA: country only","log PR ratio","Upper bound of the 95% CI from predict.rma.mv (log scale).",                                               NA_character_,
  "remra_cty_conv_ok",        "RE-MRA: country only","logical",     "TRUE if the RE-MRA model converged successfully; FALSE if convergence failed.",                              "Checked via: nlminb code == 0 (or NULL, which some metafor builds return for success), no NA or infinite beta coefficients. If FALSE, the estimates in this row are unreliable and should be treated as missing.",
  "remra_cty_status",         "RE-MRA: country only","categorical", "Outcome of the RE-MRA model fit for this country.",                                                          "ok = successful; no_conv = model did not converge (remra_cty_conv_ok = FALSE); prediction_error = model converged but predict.rma.mv failed; error = model fitting threw an error.",
  "remra_cty_prratio",        "RE-MRA: country only","PR ratio",    "RE-MRA country estimate on the PR ratio scale: exp(remra_cty_est).",                                        NA_character_,
  "remra_cty_prratio_ci_low", "RE-MRA: country only","PR ratio",    "Lower bound of the 95% CI for the RE-MRA PR ratio: exp(remra_cty_ci_low).",                                NA_character_,
  "remra_cty_prratio_ci_high","RE-MRA: country only","PR ratio",    "Upper bound of the 95% CI for the RE-MRA PR ratio: exp(remra_cty_ci_high).",                               NA_character_,

  "remra_cov_est",            "RE-MRA: covariate-adjusted","log PR ratio","RE-MRA estimate with country + eligible study-level covariates (log PR ratio). Represents expected discrimination in this country if all covariates were at their IVW-weighted sample means.", "Only fitted when degrees of freedom allow: 1 + country_df + cov_df <= n_studies (study-only DF method). Currently available for ero only; NA for gen/hed/seo/soc (df_exceeded). Covariates are scope-specific (see cty_reg_candidates_for_scope in the script).",
  "remra_cov_se",             "RE-MRA: covariate-adjusted","log PR ratio","Standard error of the covariate-adjusted RE-MRA estimate (log scale).",                              NA_character_,
  "remra_cov_ci_low",         "RE-MRA: covariate-adjusted","log PR ratio","Lower bound of the 95% CI (log scale).",                                                              NA_character_,
  "remra_cov_ci_high",        "RE-MRA: covariate-adjusted","log PR ratio","Upper bound of the 95% CI (log scale).",                                                              NA_character_,
  "remra_cov_conv_ok",        "RE-MRA: covariate-adjusted","logical",    "Convergence flag for the covariate-adjusted RE-MRA model (same logic as remra_cty_conv_ok).",          "NA when remra_cov_status = 'df_exceeded' (model was not attempted).",
  "remra_cov_status",         "RE-MRA: covariate-adjusted","categorical","Outcome of the covariate-adjusted RE-MRA fit.",                                                        "ok = successful; df_exceeded = model not attempted due to insufficient degrees of freedom; error = fitting error; no_conv = convergence failure; prediction_error = prediction step failed.",
  "remra_cov_prratio",        "RE-MRA: covariate-adjusted","PR ratio",   "Covariate-adjusted RE-MRA estimate on the PR ratio scale: exp(remra_cov_est).",                       NA_character_,
  "remra_cov_prratio_ci_low", "RE-MRA: covariate-adjusted","PR ratio",   "Lower bound of the 95% CI: exp(remra_cov_ci_low).",                                                   NA_character_,
  "remra_cov_prratio_ci_high","RE-MRA: covariate-adjusted","PR ratio",   "Upper bound of the 95% CI: exp(remra_cov_ci_high).",                                                  NA_character_,

  "uwls_cty_est",             "UWLS: country only (bias-corrected)","log PR ratio","Publication-bias-corrected country estimate from UWLS meta-regression (log PR ratio). This is the recommended primary estimate for maps.", "Estimated via fixest::feols with PET or PEESE precision term (see uwls_cty_model_type). treatment_group is absorbed as a fixed effect. Predictions use counterfactual marginal standardisation (marginaleffects::avg_predictions with grid_type = 'counterfactual'): the estimate answers 'what would average discrimination be in country X if the study-design distribution were the same as in the full dataset?' The precision term is set to SE = 0.1 at prediction time to remove the publication-bias component.",
  "uwls_cty_se",              "UWLS: country only (bias-corrected)","log PR ratio","Standard error of the UWLS country estimate (log scale).",                                  "Computed via the delta method in marginaleffects::avg_predictions.",
  "uwls_cty_ci_low",          "UWLS: country only (bias-corrected)","log PR ratio","Lower bound of the 95% CI (log scale).",                                                     NA_character_,
  "uwls_cty_ci_high",         "UWLS: country only (bias-corrected)","log PR ratio","Upper bound of the 95% CI (log scale).",                                                     NA_character_,
  "uwls_cty_conv_ok",         "UWLS: country only (bias-corrected)","logical",     "TRUE if the prediction returned a finite estimate and finite standard error.",               "FALSE indicates a numerical problem in avg_predictions; treat estimates as missing.",
  "uwls_cty_status",          "UWLS: country only (bias-corrected)","categorical", "Outcome of the UWLS country-level prediction.",                                              "ok = successful; collinear_dropped = country dummy was dropped from the model due to data sparsity — prediction would be uninformative (NA returned instead); prediction_error = avg_predictions failed; df_exceeded = model not fitted due to DF constraint; insufficient_n or insufficient_random_levels = too few observations or studies; error = model fitting error.",
  "uwls_cty_model_type",      "UWLS: country only (bias-corrected)","categorical", "Which publication-bias correction was selected for this discrimination ground.",             "pet = PET (Precision Effect Test): linear regression on SE. Selected when the PET coefficient p-value >= 0.10, indicating no statistically significant publication bias. peese = PEESE (Precision Effect Estimate with Standard Error): linear regression on variance. Selected when PET p < 0.10, indicating significant publication bias. The same model type applies to all countries within a ground. Selection mirrors the rule used in the main regression stage (Stage 5).",
  "uwls_cty_prratio",         "UWLS: country only (bias-corrected)","PR ratio",    "Bias-corrected country estimate on the PR ratio scale: exp(uwls_cty_est).",                 NA_character_,
  "uwls_cty_prratio_ci_low",  "UWLS: country only (bias-corrected)","PR ratio",    "Lower bound of the 95% CI: exp(uwls_cty_ci_low).",                                          NA_character_,
  "uwls_cty_prratio_ci_high", "UWLS: country only (bias-corrected)","PR ratio",    "Upper bound of the 95% CI: exp(uwls_cty_ci_high).",                                         NA_character_,

  "uwls_cov_est",             "UWLS: covariate-adjusted (bias-corrected)","log PR ratio","Publication-bias-corrected and covariate-adjusted country estimate (log PR ratio). Where available, this is the recommended primary estimate (preferred over uwls_cty_est).", "Same method as uwls_cty but with study-level covariates added to the model alongside country. Only fitted when DF allows. Currently available for ero only; df_exceeded for gen/hed/seo/soc. Covariates are scope-specific.",
  "uwls_cov_se",              "UWLS: covariate-adjusted (bias-corrected)","log PR ratio","Standard error (log scale).",                                                          NA_character_,
  "uwls_cov_ci_low",          "UWLS: covariate-adjusted (bias-corrected)","log PR ratio","Lower bound of the 95% CI (log scale).",                                               NA_character_,
  "uwls_cov_ci_high",         "UWLS: covariate-adjusted (bias-corrected)","log PR ratio","Upper bound of the 95% CI (log scale).",                                               NA_character_,
  "uwls_cov_conv_ok",         "UWLS: covariate-adjusted (bias-corrected)","logical",    "TRUE if the prediction returned a finite estimate and finite standard error.",           "NA when uwls_cov_status = 'df_exceeded' (model was not attempted).",
  "uwls_cov_status",          "UWLS: covariate-adjusted (bias-corrected)","categorical","Outcome of the covariate-adjusted UWLS prediction.",                                    "ok = successful; df_exceeded = model not fitted due to DF constraint (most common for gen/hed/seo/soc); none = both PET and PEESE models failed so no prediction is possible; other values same as uwls_cty_status.",
  "uwls_cov_model_type",      "UWLS: covariate-adjusted (bias-corrected)","categorical","Which bias correction was selected for the covariate-adjusted UWLS model.",             "pet / peese: same logic as uwls_cty_model_type. none = both PET and PEESE failed (typically when uwls_cov_status = 'df_exceeded' for the whole ground).",
  "uwls_cov_prratio",         "UWLS: covariate-adjusted (bias-corrected)","PR ratio",   "Bias-corrected and covariate-adjusted estimate on the PR ratio scale: exp(uwls_cov_est).", NA_character_,
  "uwls_cov_prratio_ci_low",  "UWLS: covariate-adjusted (bias-corrected)","PR ratio",   "Lower bound of the 95% CI: exp(uwls_cov_ci_low).",                                    NA_character_,
  "uwls_cov_prratio_ci_high", "UWLS: covariate-adjusted (bias-corrected)","PR ratio",   "Upper bound of the 95% CI: exp(uwls_cov_ci_high).",                                   NA_character_
)

# Write outputs ####
rds_path  <- file.path(dir.cty.out, "housing_meta_country_estimates.rds")
xlsx_path <- file.path(dir.cty.out, "housing_meta_country_estimates.xlsx")

saveRDS(housing.meta.cty.estimates, file = rds_path)

ground.scopes <- c("ero", "gen", "hed", "seo", "soc")
xlsx_sheets <- c(
  list(all = housing.meta.cty.estimates),
  purrr::set_names(
    purrr::map(ground.scopes, function(g) {
      housing.meta.cty.estimates %>% dplyr::filter(ground_abbr == g)
    }),
    ground.scopes
  ),
  list(dictionary = housing.meta.cty.dict)
)

writexl::write_xlsx(xlsx_sheets, path = xlsx_path)

# Manuscript tables ####
format_prratio_ci <- function(est, low, high) {
  ifelse(
    is.na(est),
    "\u2014",
    sprintf("%.2f [%.2f, %.2f]", est, low, high)
  )
}


ground.labels <- c(
  ero = "Ethno-racial origin",
  gen = "Gender",
  hed = "Health and disability",
  seo = "Sexual orientation",
  soc = "Social origin"
)

## Table A-C: All grounds — uwls_cty (primary) + remra_cty; k >= 2 only
tbl_all <- housing.meta.cty.estimates %>%
  dplyr::filter(k >= 2) %>%
  dplyr::filter(
    (uwls_cty_status == "ok"  & uwls_cty_conv_ok  %in% TRUE) |
    (remra_cty_status == "ok" & remra_cty_conv_ok %in% TRUE)
  ) %>%
  dplyr::mutate(
    Ground = ground.labels[ground_abbr],
    uwls_cty_fmt = format_prratio_ci(
      dplyr::if_else(uwls_cty_status == "ok" & uwls_cty_conv_ok %in% TRUE,
                     uwls_cty_prratio,         NA_real_),
      dplyr::if_else(uwls_cty_status == "ok" & uwls_cty_conv_ok %in% TRUE,
                     uwls_cty_prratio_ci_low,  NA_real_),
      dplyr::if_else(uwls_cty_status == "ok" & uwls_cty_conv_ok %in% TRUE,
                     uwls_cty_prratio_ci_high, NA_real_)
    ),
    remra_cty_fmt = format_prratio_ci(
      dplyr::if_else(remra_cty_status == "ok" & remra_cty_conv_ok %in% TRUE,
                     remra_cty_prratio,         NA_real_),
      dplyr::if_else(remra_cty_status == "ok" & remra_cty_conv_ok %in% TRUE,
                     remra_cty_prratio_ci_low,  NA_real_),
      dplyr::if_else(remra_cty_status == "ok" & remra_cty_conv_ok %in% TRUE,
                     remra_cty_prratio_ci_high, NA_real_)
    )
  ) %>%
  dplyr::arrange(Ground, country) %>%
  dplyr::select(
    Ground,
    Country                          = country,
    k                                = k,
    Studies                          = n_studies,
    `UWLS (bias-corrected) [95% CI]` = uwls_cty_fmt,
    `RE-MRA [95% CI]`                = remra_cty_fmt
  )

xlsx_tables_path <- file.path(dir.cty.out, "meta_housing_country_estimates_tables.xlsx")
writexl::write_xlsx(list(`Table A-C` = tbl_all), path = xlsx_tables_path)

message("Country estimates complete.")
