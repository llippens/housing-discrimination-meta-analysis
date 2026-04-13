# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "dplyr", "purrr", "tibble",
  "metafor", "clubSandwich", "stringr", "writexl"
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
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_specs.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_output.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reporting_workbooks.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.meta.reg <- file.path(dir.tables, "3_reg_remra")

ensure_dir(dir.tables)
dir.meta.reg.dirs <- reg_output_dirs(dir.meta.reg)
dir.meta.reg.appendix <- dir.meta.reg.dirs$appendix
dir.meta.reg.appendix.sensitivity <- dir.meta.reg.dirs$sensitivity
dir.meta.reg.appendix.interactions <- dir.meta.reg.dirs$interactions

# Cleanup ####
stale.paths <- c(
  Sys.glob(file.path(dir.meta.reg, "meta_housing_reg_ground_*.csv")),
  Sys.glob(file.path(dir.meta.reg, "housing_meta_reg_ground_*.rds")),
  file.path(dir.meta.reg, "meta_housing_reg_by_ground.csv"),
  file.path(dir.meta.reg, "housing_meta_reg_by_ground_list.rds")
)
stale.paths <- stale.paths[file.exists(stale.paths)]
if (length(stale.paths) > 0) {
  file.remove(stale.paths)
}

stale.sensitivity.paths <- c(
  Sys.glob(file.path(dir.meta.reg.appendix.sensitivity, "meta_housing_reg_*_ymf*.csv")),
  Sys.glob(file.path(dir.meta.reg.appendix.sensitivity, "housing_meta_reg_*_ymf*.rds")),
  Sys.glob(file.path(dir.meta.reg.appendix.sensitivity, "meta_housing_reg_*_ctr*.csv")),
  Sys.glob(file.path(dir.meta.reg.appendix.sensitivity, "housing_meta_reg_*_ctr*.rds")),
  file.path(dir.meta.reg.appendix.sensitivity, "housing_meta_reg_ymf.rds"),
  file.path(dir.meta.reg.appendix.sensitivity, "housing_meta_reg_ctr.rds")
 )
stale.sensitivity.paths <- stale.sensitivity.paths[file.exists(stale.sensitivity.paths)]
if (length(stale.sensitivity.paths) > 0) {
  file.remove(stale.sensitivity.paths)
}

stale.interaction.paths <- c(
  Sys.glob(file.path(dir.meta.reg.appendix.interactions, "meta_housing_reg_*_int_*.csv")),
  Sys.glob(file.path(dir.meta.reg.appendix.interactions, "housing_meta_reg_*_int_*.rds"))
)
stale.interaction.paths <- stale.interaction.paths[file.exists(stale.interaction.paths)]
if (length(stale.interaction.paths) > 0) {
  file.remove(stale.interaction.paths)
}

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/0_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "prratio_ln", "prratio_ln_v", "study", "cluster", "ground_abbr",
  "region_agg", "country", "year_mid", "year_mid_fct", "info_binary", "agent_type",
  "matched_design", "audit_type", "treatment_group", "callback_type",
  "airbnb", "can_gender", "can_education", "can_employment",
  "can_migrant_generation", "peer_reviewed", "language"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Model setup ####
remra.specs <- reg_specs_remra()
factor.cols <- c(
  "region_agg", "country", "info_binary", "agent_type", "matched_design", "audit_type",
  "treatment_group", "callback_type", "airbnb", "can_gender",
  "can_education", "can_employment", "can_migrant_generation",
  "peer_reviewed", "language", "period_narrow"
)
numeric.cols <- c("year_mid")
factor.cols.sens <- unique(c(factor.cols, "year_mid_fct"))

# Model runs ####
housing.meta.by.ground <- split(housing.meta, as.character(housing.meta$ground_abbr))
housing.meta.scope.data <- c(list(overall = housing.meta), housing.meta.by.ground)

housing.meta.design.list <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    build_metareg_design(
      data = df,
      scope_name = scope_name,
      factor_cols = factor.cols,
      numeric_cols = numeric.cols
    )
  }
)
housing.meta.design <- purrr::map_dfr(housing.meta.design.list, "design")

housing.meta.fit.list <- purrr::imap(
  housing.meta.design.list,
  function(scope_obj, scope_name) {
    fit_metareg_scope(scope_obj, scope_name, factor.cols)
  }
)

housing.meta.table.overall <- housing.meta.fit.list$overall$table %>%
  mutate(model_scope = "overall")

housing.meta.tables.by.ground <- housing.meta.fit.list[names(housing.meta.fit.list) != "overall"] %>%
  purrr::imap(
    function(obj, g) {
      obj$table %>%
        mutate(
          model_scope = paste0("ground_", g),
          ground_abbr = g
        )
    }
  )

housing.meta.fit.status <- purrr::map_dfr(housing.meta.fit.list, "fit_status") %>%
  mutate(scope = as.character(scope))

housing.meta.fit.status <- add_pseudo_r2_to_fit_status(
  housing.meta.design.list,
  housing.meta.fit.list,
  housing.meta.fit.status
)

# Outputs ####
housing.meta.reg <- list(
  housing.meta.design = housing.meta.design,
  housing.meta.fit.status = housing.meta.fit.status,
  housing.meta.model.overall = housing.meta.fit.list$overall,
  housing.meta.models.by.ground = housing.meta.fit.list[names(housing.meta.fit.list) != "overall"],
  housing.meta.table.overall = housing.meta.table.overall,
  housing.meta.tables.by.ground = housing.meta.tables.by.ground
)

write_remra_output_bundle(
  reg_obj = housing.meta.reg,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$main),
  tag = remra.specs$main$tag
)

# Sensitivity: Year midpoint factor ####
housing.meta.design.list.ymf <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    fixed_candidates_sens <- fixed_candidates_for_scope(scope_name)
    fixed_candidates_sens[fixed_candidates_sens == "year_mid"] <- "year_mid_fct"
    fixed_candidates_sens <- unique(fixed_candidates_sens)

    build_metareg_design(
      data = df,
      scope_name = scope_name,
      factor_cols = factor.cols.sens,
      numeric_cols = numeric.cols,
      fixed_candidates_override = fixed_candidates_sens
    )
  }
)
housing.meta.design.ymf <- purrr::map_dfr(housing.meta.design.list.ymf, "design")

housing.meta.fit.list.ymf <- purrr::imap(
  housing.meta.design.list.ymf,
  function(scope_obj, scope_name) {
    fit_metareg_scope(scope_obj, scope_name, factor.cols.sens)
  }
)

housing.meta.table.overall.ymf <- housing.meta.fit.list.ymf$overall$table %>%
  mutate(model_scope = "overall")

housing.meta.tables.by.ground.ymf <- housing.meta.fit.list.ymf[names(housing.meta.fit.list.ymf) != "overall"] %>%
  purrr::imap(
    function(obj, g) {
      obj$table %>%
        mutate(
          model_scope = paste0("ground_", g),
          ground_abbr = g
        )
    }
  )

housing.meta.fit.status.ymf <- purrr::map_dfr(housing.meta.fit.list.ymf, "fit_status") %>%
  mutate(scope = as.character(scope))

housing.meta.fit.status.ymf <- add_pseudo_r2_to_fit_status(
  housing.meta.design.list.ymf,
  housing.meta.fit.list.ymf,
  housing.meta.fit.status.ymf
)

housing.meta.reg.ymf <- list(
  housing.meta.design = housing.meta.design.ymf,
  housing.meta.fit.status = housing.meta.fit.status.ymf,
  housing.meta.model.overall = housing.meta.fit.list.ymf$overall,
  housing.meta.models.by.ground = housing.meta.fit.list.ymf[names(housing.meta.fit.list.ymf) != "overall"],
  housing.meta.table.overall = housing.meta.table.overall.ymf,
  housing.meta.tables.by.ground = housing.meta.tables.by.ground.ymf
)

write_remra_output_bundle(
  reg_obj = housing.meta.reg.ymf,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$ymf),
  tag = remra.specs$ymf$tag
)

# Sensitivity: Period (5-year blocks) ####
housing.meta.design.list.period <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    fixed_candidates_period <- fixed_candidates_for_scope(scope_name)
    fixed_candidates_period[fixed_candidates_period == "year_mid"] <- "period_narrow"
    fixed_candidates_period <- unique(fixed_candidates_period)

    build_metareg_design(
      data = df,
      scope_name = scope_name,
      factor_cols = factor.cols,
      numeric_cols = numeric.cols,
      fixed_candidates_override = fixed_candidates_period
    )
  }
)
housing.meta.design.period <- purrr::map_dfr(housing.meta.design.list.period, "design")

housing.meta.fit.list.period <- purrr::imap(
  housing.meta.design.list.period,
  function(scope_obj, scope_name) {
    fit_metareg_scope(scope_obj, scope_name, factor.cols)
  }
)

housing.meta.table.overall.period <- housing.meta.fit.list.period$overall$table %>%
  mutate(model_scope = "overall")

housing.meta.tables.by.ground.period <- housing.meta.fit.list.period[names(housing.meta.fit.list.period) != "overall"] %>%
  purrr::imap(function(obj, g) {
    obj$table %>% mutate(model_scope = paste0("ground_", g), ground_abbr = g)
  })

housing.meta.fit.status.period <- purrr::map_dfr(housing.meta.fit.list.period, "fit_status") %>%
  mutate(scope = as.character(scope))

housing.meta.fit.status.period <- add_pseudo_r2_to_fit_status(
  housing.meta.design.list.period,
  housing.meta.fit.list.period,
  housing.meta.fit.status.period
)

housing.meta.reg.period <- list(
  housing.meta.design           = housing.meta.design.period,
  housing.meta.fit.status       = housing.meta.fit.status.period,
  housing.meta.model.overall    = housing.meta.fit.list.period$overall,
  housing.meta.models.by.ground = housing.meta.fit.list.period[names(housing.meta.fit.list.period) != "overall"],
  housing.meta.table.overall    = housing.meta.table.overall.period,
  housing.meta.tables.by.ground = housing.meta.tables.by.ground.period
)

write_remra_output_bundle(
  reg_obj    = housing.meta.reg.period,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$period),
  tag        = remra.specs$period$tag
)

# Sensitivity: Country regressor ####
housing.meta.design.list.ctr <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    fixed_candidates_sens <- fixed_candidates_for_scope(scope_name)
    fixed_candidates_sens[fixed_candidates_sens == "region_agg"] <- "country"
    if (!"country" %in% fixed_candidates_sens) {
      fixed_candidates_sens <- c(fixed_candidates_sens, "country")
    }
    fixed_candidates_sens <- unique(fixed_candidates_sens)

    build_metareg_design(
      data = df,
      scope_name = scope_name,
      factor_cols = factor.cols,
      numeric_cols = numeric.cols,
      fixed_candidates_override = fixed_candidates_sens
    )
  }
)
housing.meta.design.ctr <- purrr::map_dfr(housing.meta.design.list.ctr, "design")

housing.meta.fit.list.ctr <- purrr::imap(
  housing.meta.design.list.ctr,
  function(scope_obj, scope_name) {
    fit_metareg_scope(scope_obj, scope_name, factor.cols)
  }
)

housing.meta.table.overall.ctr <- housing.meta.fit.list.ctr$overall$table %>%
  mutate(model_scope = "overall")

housing.meta.tables.by.ground.ctr <- housing.meta.fit.list.ctr[names(housing.meta.fit.list.ctr) != "overall"] %>%
  purrr::imap(
    function(obj, g) {
      obj$table %>%
        mutate(
          model_scope = paste0("ground_", g),
          ground_abbr = g
        )
    }
  )

housing.meta.fit.status.ctr <- purrr::map_dfr(housing.meta.fit.list.ctr, "fit_status") %>%
  mutate(scope = as.character(scope))

housing.meta.fit.status.ctr <- add_pseudo_r2_to_fit_status(
  housing.meta.design.list.ctr,
  housing.meta.fit.list.ctr,
  housing.meta.fit.status.ctr
)

housing.meta.reg.ctr <- list(
  housing.meta.design = housing.meta.design.ctr,
  housing.meta.fit.status = housing.meta.fit.status.ctr,
  housing.meta.model.overall = housing.meta.fit.list.ctr$overall,
  housing.meta.models.by.ground = housing.meta.fit.list.ctr[names(housing.meta.fit.list.ctr) != "overall"],
  housing.meta.table.overall = housing.meta.table.overall.ctr,
  housing.meta.tables.by.ground = housing.meta.tables.by.ground.ctr
)

write_remra_output_bundle(
  reg_obj = housing.meta.reg.ctr,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$ctr),
  tag = remra.specs$ctr$tag
)

run_metareg_sensitivity <- function(
    tag,
    output_dir,
    factor_cols_use,
    numeric_cols_use,
    fixed_candidates_fn = NULL,
    interaction_terms_fn = NULL
) {
  design_list <- purrr::imap(
    housing.meta.scope.data,
    function(df, scope_name) {
      fixed_candidates <- if (is.null(fixed_candidates_fn)) {
        fixed_candidates_for_scope(scope_name)
      } else {
        fixed_candidates_fn(scope_name)
      }
      interaction_terms <- if (is.null(interaction_terms_fn)) {
        character(0)
      } else {
        interaction_terms_fn(scope_name)
      }

      build_metareg_design(
        data = df,
        scope_name = scope_name,
        factor_cols = factor_cols_use,
        numeric_cols = numeric_cols_use,
        fixed_candidates_override = fixed_candidates,
        interaction_terms = interaction_terms
      )
    }
  )
  design_tbl <- purrr::map_dfr(design_list, "design")

  fit_list <- purrr::imap(
    design_list,
    function(scope_obj, scope_name) {
      fit_metareg_scope(scope_obj, scope_name, factor_cols_use)
    }
  )

  table_overall <- fit_list$overall$table %>%
    mutate(model_scope = "overall")

  tables_by_ground <- fit_list[names(fit_list) != "overall"] %>%
    purrr::imap(
      function(obj, g) {
        obj$table %>%
          mutate(
            model_scope = paste0("ground_", g),
            ground_abbr = g
          )
      }
    )

  fit_status_tbl <- purrr::map_dfr(fit_list, "fit_status") %>%
    mutate(scope = as.character(scope))

  fit_status_tbl <- add_pseudo_r2_to_fit_status(design_list, fit_list, fit_status_tbl)

  reg_obj <- list(
    housing.meta.design = design_tbl,
    housing.meta.fit.status = fit_status_tbl,
    housing.meta.model.overall = fit_list$overall,
    housing.meta.models.by.ground = fit_list[names(fit_list) != "overall"],
    housing.meta.table.overall = table_overall,
    housing.meta.tables.by.ground = tables_by_ground
  )

  write_remra_output_bundle(
    reg_obj = reg_obj,
    output_dir = output_dir,
    tag = tag
  )

  return(reg_obj)
}

# Sensitivity: Region x treatment group interaction ####
housing.meta.reg.int.region.tg <- run_metareg_sensitivity(
  tag = remra.specs$int_region_tg$tag,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$int_region_tg),
  factor_cols_use = factor.cols,
  numeric_cols_use = numeric.cols,
  fixed_candidates_fn = function(scope_name) fixed_candidates_for_scope(scope_name),
  interaction_terms_fn = function(scope_name) "region_agg:treatment_group"
)

# Sensitivity: Region x year interaction ####
housing.meta.reg.int.region.year <- run_metareg_sensitivity(
  tag = remra.specs$int_region_year$tag,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$int_region_year),
  factor_cols_use = factor.cols,
  numeric_cols_use = numeric.cols,
  fixed_candidates_fn = function(scope_name) fixed_candidates_for_scope(scope_name),
  interaction_terms_fn = function(scope_name) "region_agg:year_mid"
)

# Sensitivity: Region x period interaction ####
housing.meta.reg.int.region.period <- run_metareg_sensitivity(
  tag = remra.specs$int_region_period$tag,
  output_dir = reg_output_dir_for_spec(dir.meta.reg.dirs, remra.specs$int_region_period),
  factor_cols_use = factor.cols,
  numeric_cols_use = numeric.cols,
  fixed_candidates_fn = function(scope_name) {
    cands <- fixed_candidates_for_scope(scope_name)
    cands[cands == "year_mid"] <- "period_narrow"
    unique(cands)
  },
  interaction_terms_fn = function(scope_name) "region_agg:period_narrow"
)

export_remra_reporting_workbooks(
  dir_map = dir.meta.reg.dirs,
  specs = remra.specs,
  reg_objects = list(
    main = housing.meta.reg,
    ymf = housing.meta.reg.ymf,
    ctr = housing.meta.reg.ctr,
    int_region_tg = housing.meta.reg.int.region.tg,
    int_region_year = housing.meta.reg.int.region.year,
    period = housing.meta.reg.period,
    int_region_period = housing.meta.reg.int.region.period
  )
)

message("Meta-regression complete.")
