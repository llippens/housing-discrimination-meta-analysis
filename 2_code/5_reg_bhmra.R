# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c("here", "dplyr", "tibble", "purrr", "writexl")
missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_design_utils.R"))
source(file.path(here::here(), "2_code", "functions", "reg", "f_bayes_helpers.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_specs.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_output.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reporting_workbooks.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.bayesian <- file.path(dir.tables, "5_reg_bhmra")

ensure_dir(dir.tables)
dir.bayesian.dirs <- reg_output_dirs(dir.bayesian)
dir.bayesian.appendix <- dir.bayesian.dirs$appendix
dir.bayesian.appendix.sensitivity <- dir.bayesian.dirs$sensitivity
dir.bayesian.appendix.interactions <- dir.bayesian.dirs$interactions

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/0_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)
assert_required_cols(
  housing.meta,
  c(
    "prratio_ln", "prratio_ln_v", "study", "country", "year_mid_fct", "treatment_group",
    "ground_abbr",
    "region_agg", "year_mid", "info_binary", "agent_type",
    "matched_design", "audit_type", "callback_type", "airbnb",
    "can_gender", "can_education", "can_employment", "can_migrant_generation",
    "peer_reviewed", "language"
  ),
  "housing.meta data"
)

# Model setup ####
bhma.specs <- reg_specs_bhmra()
random.vars.main <- c("study", "treatment_group")
random.vars.ctr <- c("study", "treatment_group")
random.vars.int.region.tg <- c("study")
random.vars.int.region.year <- c("study", "treatment_group")
factor.cols <- c(
  "study", "country", "year_mid_fct", "treatment_group",
  "ground_abbr",
  "region_agg", "info_binary", "agent_type",
  "matched_design", "audit_type", "callback_type", "airbnb",
  "can_gender", "can_education", "can_employment", "can_migrant_generation",
  "peer_reviewed", "language", "period_narrow"
)

# Run controls ####
parse_bool_flag <- function(value, default = FALSE) {
  if (is.null(value) || length(value) == 0 || is.na(value) || value == "") {
    return(default)
  }
  value_norm <- tolower(trimws(as.character(value[[1]])))
  if (value_norm %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value_norm %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  default
}

get_cli_flag <- function(args, name, default = FALSE) {
  pattern <- paste0("^--", name, "=")
  matches <- args[grepl(pattern, args)]
  if (length(matches) == 0) {
    return(default)
  }
  value <- sub(pattern, "", matches[[length(matches)]])
  parse_bool_flag(value, default = default)
}

cli.args <- commandArgs(trailingOnly = TRUE)
run.main <- TRUE
run.optional.all <- get_cli_flag(cli.args, "run_optional_all", default = FALSE)
run.ctr <- get_cli_flag(cli.args, "run_ctr", default = FALSE) || run.optional.all
run.int.region.tg <- get_cli_flag(cli.args, "run_int_region_tg", default = FALSE) || run.optional.all
run.int.region.year <- get_cli_flag(cli.args, "run_int_region_year", default = FALSE) || run.optional.all
run.sens.period       <- get_cli_flag(cli.args, "run_sens_period",       default = FALSE) || run.optional.all
run.int.region.period <- get_cli_flag(cli.args, "run_int_region_period", default = FALSE) || run.optional.all

run.config <- tibble(
  run_main = run.main,
  run_ctr = run.ctr,
  run_int_region_tg = run.int.region.tg,
  run_int_region_year = run.int.region.year,
  run_sens_period = run.sens.period,
  run_int_region_period = run.int.region.period,
  run_optional_all = run.optional.all,
  cli_args = if (length(cli.args) == 0) NA_character_ else paste(cli.args, collapse = " ")
)
saveRDS(run.config, file = file.path(dir.bayesian, "housing_meta_bayes_run_config.rds"))
message(
  "BHM run controls -> main: ", run.main,
  ", ctr: ", run.ctr,
  ", int_region_tg: ", run.int.region.tg,
  ", int_region_year: ", run.int.region.year,
  ", sens_period: ", run.sens.period,
  ", int_region_period: ", run.int.region.period,
  "."
)

# Helper functions ####
fixed_terms_for_scope <- function(
    scope_name,
    geography_term = "region_agg",
    include_treatment_group = FALSE
) {
  terms <- c(
    geography_term, "year_mid", "info_binary", "agent_type",
    "matched_design", "audit_type", "callback_type", "airbnb",
    "can_education", "can_employment", "peer_reviewed", "language"
  )
  if (include_treatment_group) {
    terms <- c(terms, "treatment_group")
  }
  if (scope_name != "gen") {
    terms <- c(terms, "can_gender")
  }
  if (scope_name == "ero") {
    terms <- c(terms, "can_migrant_generation")
  }

  unique(terms)
}

# Design ####
housing.meta.by.ground <- split(housing.meta, as.character(housing.meta$ground_abbr))
housing.meta.model.sets <- c(list(overall = housing.meta), housing.meta.by.ground)

housing.meta.bayes.design.list <- purrr::imap(
  housing.meta.model.sets,
  function(df, scope_name) {
    build_bayes_design(
      data = df,
      scope_name = scope_name,
      fixed_candidates = fixed_terms_for_scope(scope_name, geography_term = "region_agg"),
      random_terms = random.vars.main,
      factor_cols = factor.cols,
      interaction_terms = character(0)
    )
  }
)
housing.meta.bayes.design <- purrr::map_dfr(housing.meta.bayes.design.list, "design")

saveRDS(housing.meta.bayes.design, file = file.path(dir.bayesian, "housing_meta_bayes_design_check.rds"))

housing.meta.bayes.design.list.ctr <- NULL
housing.meta.bayes.design.ctr <- NULL
if (run.ctr) {
  housing.meta.bayes.design.list.ctr <- purrr::imap(
    housing.meta.model.sets,
    function(df, scope_name) {
      build_bayes_design(
        data = df,
        scope_name = scope_name,
        fixed_candidates = fixed_terms_for_scope(scope_name, geography_term = "country"),
        random_terms = random.vars.ctr,
        factor_cols = factor.cols,
        interaction_terms = character(0)
      )
    }
  )
  housing.meta.bayes.design.ctr <- purrr::map_dfr(housing.meta.bayes.design.list.ctr, "design")
  saveRDS(
    housing.meta.bayes.design.ctr,
    file = file.path(
      reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$ctr),
      paste0("housing_meta_bayes_design_check_", bhma.specs$ctr$tag, ".rds")
    )
  )
} else {
  message("Skipping country sensitivity design (run_ctr = FALSE).")
}

housing.meta.bayes.design.list.int.region.tg <- NULL
housing.meta.bayes.design.int.region.tg <- NULL
if (run.int.region.tg) {
  housing.meta.bayes.design.list.int.region.tg <- purrr::imap(
    housing.meta.model.sets,
    function(df, scope_name) {
      build_bayes_design(
        data = df,
        scope_name = scope_name,
        fixed_candidates = fixed_terms_for_scope(
          scope_name,
          geography_term = "region_agg",
          include_treatment_group = TRUE
        ),
        random_terms = random.vars.int.region.tg,
        factor_cols = factor.cols,
        interaction_terms = c("region_agg:treatment_group")
      )
    }
  )
  housing.meta.bayes.design.int.region.tg <- purrr::map_dfr(
    housing.meta.bayes.design.list.int.region.tg,
    "design"
  )
  saveRDS(
    housing.meta.bayes.design.int.region.tg,
    file = file.path(
      reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_tg),
      paste0("housing_meta_bayes_design_check_", bhma.specs$int_region_tg$tag, ".rds")
    )
  )
} else {
  message("Skipping region x treatment-group interaction design (run_int_region_tg = FALSE).")
}

housing.meta.bayes.design.list.int.region.year <- NULL
housing.meta.bayes.design.int.region.year <- NULL
if (run.int.region.year) {
  housing.meta.bayes.design.list.int.region.year <- purrr::imap(
    housing.meta.model.sets,
    function(df, scope_name) {
      build_bayes_design(
        data = df,
        scope_name = scope_name,
        fixed_candidates = fixed_terms_for_scope(scope_name, geography_term = "region_agg"),
        random_terms = random.vars.int.region.year,
        factor_cols = factor.cols,
        interaction_terms = c("region_agg:year_mid")
      )
    }
  )
  housing.meta.bayes.design.int.region.year <- purrr::map_dfr(
    housing.meta.bayes.design.list.int.region.year,
    "design"
  )
  saveRDS(
    housing.meta.bayes.design.int.region.year,
    file = file.path(
      reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_year),
      paste0("housing_meta_bayes_design_check_", bhma.specs$int_region_year$tag, ".rds")
    )
  )
} else {
  message("Skipping region x year interaction design (run_int_region_year = FALSE).")
}

housing.meta.bayes.design.list.period <- NULL
housing.meta.bayes.design.period <- NULL
if (run.sens.period) {
  housing.meta.bayes.design.list.period <- purrr::imap(
    housing.meta.model.sets,
    function(df, scope_name) {
      fixed_candidates_period <- fixed_terms_for_scope(scope_name, geography_term = "region_agg")
      fixed_candidates_period[fixed_candidates_period == "year_mid"] <- "period_narrow"
      fixed_candidates_period <- unique(fixed_candidates_period)

      build_bayes_design(
        data = df,
        scope_name = scope_name,
        fixed_candidates = fixed_candidates_period,
        random_terms = random.vars.main,
        factor_cols = factor.cols,
        interaction_terms = character(0)
      )
    }
  )
  housing.meta.bayes.design.period <- purrr::map_dfr(
    housing.meta.bayes.design.list.period,
    "design"
  )
  saveRDS(
    housing.meta.bayes.design.period,
    file = file.path(
      reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$period),
      paste0("housing_meta_bayes_design_check_", bhma.specs$period$tag, ".rds")
    )
  )
} else {
  message("Skipping period sensitivity design (run_sens_period = FALSE).")
}

housing.meta.bayes.design.list.int.region.period <- NULL
housing.meta.bayes.design.int.region.period <- NULL
if (run.int.region.period) {
  housing.meta.bayes.design.list.int.region.period <- purrr::imap(
    housing.meta.model.sets,
    function(df, scope_name) {
      fixed_candidates_period <- fixed_terms_for_scope(scope_name, geography_term = "region_agg")
      fixed_candidates_period[fixed_candidates_period == "year_mid"] <- "period_narrow"
      fixed_candidates_period <- unique(fixed_candidates_period)

      build_bayes_design(
        data = df,
        scope_name = scope_name,
        fixed_candidates = fixed_candidates_period,
        random_terms = random.vars.int.region.year,
        factor_cols = factor.cols,
        interaction_terms = c("region_agg:period_narrow")
      )
    }
  )
  housing.meta.bayes.design.int.region.period <- purrr::map_dfr(
    housing.meta.bayes.design.list.int.region.period,
    "design"
  )
  saveRDS(
    housing.meta.bayes.design.int.region.period,
    file = file.path(
      reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_period),
      paste0("housing_meta_bayes_design_check_", bhma.specs$int_region_period$tag, ".rds")
    )
  )
} else {
  message("Skipping region x period interaction design (run_int_region_period = FALSE).")
}

# Bayesian fits ####
if (!requireNamespace("brms", quietly = TRUE)) {
  stop("Package brms is required for 2_code/5_reg_bhmra.R.", call. = FALSE)
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Package cmdstanr is required for 2_code/5_reg_bhmra.R.", call. = FALSE)
}

bayes.backend <- "cmdstanr"
priors <- c(
  brms::prior(normal(0, 1), class = "Intercept"),
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(cauchy(0, 0.5), class = "sd"),
  brms::prior(cauchy(0, 0.5), class = "sigma")
)

fit_bayes_set <- function(design_list, priors, bayes.backend) {
  models <- list()
  fit_status <- list()
  summary <- list()

  for (scope_name in names(design_list)) {
    scope_obj <- design_list[[scope_name]]
    if (scope_obj$status != "ready") {
      fit_status[[scope_name]] <- tibble(
        scope = scope_name,
        fit_status = "skipped",
        fit_note = scope_obj$status
      )
      next
    }

    scope_formula <- build_bayes_formula(
      fixed_terms = scope_obj$fixed_terms,
      interaction_terms = scope_obj$interaction_terms,
      random_terms = scope_obj$random_terms
    )
    scope_fit <- tryCatch(
      brms::brm(
        formula = scope_formula,
        data = scope_obj$data,
        family = brms::student(),
        prior = priors,
        iter = 10000,
        warmup = 2000,
        chains = 4,
        cores = min(4, parallel::detectCores(logical = FALSE)),
        seed = 8888,
        backend = bayes.backend,
        refresh = 0,
        control = list(adapt_delta = 0.95, max_treedepth = 15)
      ),
      error = function(e) e
    )

    if (inherits(scope_fit, "error")) {
      fit_status[[scope_name]] <- tibble(
        scope = scope_name,
        fit_status = "error",
        fit_note = conditionMessage(scope_fit)
      )
      next
    }

    models[[scope_name]] <- scope_fit
    fit_status[[scope_name]] <- tibble(
      scope = scope_name,
      fit_status = "ok",
      fit_note = NA_character_
    )
    summary[[scope_name]] <- extract_bayes_summary(scope_fit, scope_name)
  }

  fit_status <- bind_rows(fit_status)

  if (length(summary) == 0) {
    summary <- tibble(
      term = character(0),
      estimate = numeric(0),
      est.error = numeric(0),
      q2.5 = numeric(0),
      q97.5 = numeric(0),
      scope = character(0),
      component = character(0),
      n = integer(0),
      j_study = integer(0),
      j_treatment_group = integer(0)
    )
  } else {
    summary <- bind_rows(summary)
  }

  status <- tibble(
    status = "completed",
    backend = bayes.backend,
    timestamp = as.character(Sys.time()),
    n_models_ok = sum(fit_status$fit_status == "ok"),
    n_models_total = nrow(fit_status)
  )

  list(
    models = models,
    fit_status = fit_status,
    summary = summary,
    status = status
  )
}

bayes.main <- fit_bayes_set(housing.meta.bayes.design.list, priors, bayes.backend)

bayes.ctr <- NULL
if (run.ctr) {
  bayes.ctr <- fit_bayes_set(housing.meta.bayes.design.list.ctr, priors, bayes.backend)
}

bayes.int.region.tg <- NULL
if (run.int.region.tg) {
  bayes.int.region.tg <- fit_bayes_set(housing.meta.bayes.design.list.int.region.tg, priors, bayes.backend)
}

bayes.int.region.year <- NULL
if (run.int.region.year) {
  bayes.int.region.year <- fit_bayes_set(housing.meta.bayes.design.list.int.region.year, priors, bayes.backend)
}

bayes.sens.period <- NULL
if (run.sens.period) {
  bayes.sens.period <- fit_bayes_set(housing.meta.bayes.design.list.period, priors, bayes.backend)
}

bayes.int.region.period <- NULL
if (run.int.region.period) {
  bayes.int.region.period <- fit_bayes_set(housing.meta.bayes.design.list.int.region.period, priors, bayes.backend)
}

# Outputs ####
write_bayes_output_bundle(
  bayes_obj = bayes.main,
  output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$main),
  tag = bhma.specs$main$tag
)

if (run.ctr) {
  write_bayes_output_bundle(
    bayes_obj = bayes.ctr,
    output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$ctr),
    tag = bhma.specs$ctr$tag
  )
}
if (run.int.region.tg) {
  write_bayes_output_bundle(
    bayes_obj = bayes.int.region.tg,
    output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_tg),
    tag = bhma.specs$int_region_tg$tag
  )
}
if (run.int.region.year) {
  write_bayes_output_bundle(
    bayes_obj = bayes.int.region.year,
    output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_year),
    tag = bhma.specs$int_region_year$tag
  )
}
if (run.sens.period) {
  write_bayes_output_bundle(
    bayes_obj  = bayes.sens.period,
    output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$period),
    tag        = bhma.specs$period$tag
  )
}
if (run.int.region.period) {
  write_bayes_output_bundle(
    bayes_obj  = bayes.int.region.period,
    output_dir = reg_output_dir_for_spec(dir.bayesian.dirs, bhma.specs$int_region_period),
    tag        = bhma.specs$int_region_period$tag
  )
}

export_bhm_reporting_workbooks(
  dir_map = dir.bayesian.dirs,
  specs = bhma.specs,
  bayes_objects = list(
    main = bayes.main,
    ctr = bayes.ctr,
    int_region_tg = bayes.int.region.tg,
    int_region_year = bayes.int.region.year,
    period = bayes.sens.period,
    int_region_period = bayes.int.region.period
  ),
  design_tables = list(
    main = housing.meta.bayes.design,
    ctr = housing.meta.bayes.design.ctr,
    int_region_tg = housing.meta.bayes.design.int.region.tg,
    int_region_year = housing.meta.bayes.design.int.region.year,
    period = housing.meta.bayes.design.period,
    int_region_period = housing.meta.bayes.design.int.region.period
  )
)

message(
  "Bayesian models completed (main): ",
  sum(bayes.main$fit_status$fit_status == "ok"),
  " / ",
  nrow(bayes.main$fit_status),
  "."
)
message(
  "Bayesian models completed (country sensitivity): ",
  if (run.ctr) {
    paste0(sum(bayes.ctr$fit_status$fit_status == "ok"), " / ", nrow(bayes.ctr$fit_status))
  } else {
    "skipped"
  },
  "."
)
message(
  "Bayesian models completed (region x treatment interaction): ",
  if (run.int.region.tg) {
    paste0(sum(bayes.int.region.tg$fit_status$fit_status == "ok"), " / ", nrow(bayes.int.region.tg$fit_status))
  } else {
    "skipped"
  },
  "."
)
message(
  "Bayesian models completed (region x year interaction): ",
  if (run.int.region.year) {
    paste0(sum(bayes.int.region.year$fit_status$fit_status == "ok"), " / ", nrow(bayes.int.region.year$fit_status))
  } else {
    "skipped"
  },
  "."
)
message(
  "Bayesian models completed (period sensitivity): ",
  if (run.sens.period) {
    paste0(sum(bayes.sens.period$fit_status$fit_status == "ok"), " / ", nrow(bayes.sens.period$fit_status))
  } else {
    "skipped"
  },
  "."
)
message(
  "Bayesian models completed (region x period interaction): ",
  if (run.int.region.period) {
    paste0(sum(bayes.int.region.period$fit_status$fit_status == "ok"), " / ", nrow(bayes.int.region.period$fit_status))
  } else {
    "skipped"
  },
  "."
)
message("BHM complete.")
