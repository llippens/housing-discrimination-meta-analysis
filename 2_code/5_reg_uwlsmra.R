# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "dplyr", "purrr", "tibble",
  "fixest", "stringr", "writexl"
)

missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_design_utils.R"))
source(file.path(here::here(), "2_code", "functions", "reg", "f_uwls_helpers.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_specs.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_output.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reporting_workbooks.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.uwls <- file.path(dir.tables, "5_reg_uwlsmra")

ensure_dir(dir.tables)
dir.uwls.dirs <- reg_output_dirs(dir.uwls)
dir.uwls.appendix <- dir.uwls.dirs$appendix
dir.uwls.appendix.sensitivity <- dir.uwls.dirs$sensitivity
dir.uwls.appendix.interactions <- dir.uwls.dirs$interactions

# Cleanup ####
stale.sensitivity.paths <- c(
  Sys.glob(file.path(dir.uwls.appendix.sensitivity, "meta_housing_uwls_*_sens_*.csv")),
  Sys.glob(file.path(dir.uwls.appendix.sensitivity, "housing_meta_uwls_models_sens_*.rds"))
)
stale.sensitivity.paths <- stale.sensitivity.paths[file.exists(stale.sensitivity.paths)]
if (length(stale.sensitivity.paths) > 0) {
  file.remove(stale.sensitivity.paths)
}

stale.interaction.paths <- c(
  Sys.glob(file.path(dir.uwls.appendix.interactions, "meta_housing_uwls_*_int_*.csv")),
  Sys.glob(file.path(dir.uwls.appendix.interactions, "housing_meta_uwls_models_int_*.rds"))
)
stale.interaction.paths <- stale.interaction.paths[file.exists(stale.interaction.paths)]
if (length(stale.interaction.paths) > 0) {
  file.remove(stale.interaction.paths)
}

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/1_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "prratio_ln", "prratio_ln_v", "study", "ground_abbr",
  "country", "treatment_group", "region_agg", "year_mid", "year_mid_fct",
  "info_binary", "agent_type", "matched_design", "audit_type",
  "callback_type", "airbnb", "can_gender", "can_education",
  "can_employment", "can_migrant_generation", "peer_reviewed", "language"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Model setup ####
uwls.specs <- reg_specs_uwls()
factor.cols <- c(
  "region_agg", "country", "info_binary", "agent_type", "matched_design", "audit_type",
  "callback_type", "airbnb", "can_gender", "can_education",
  "can_employment", "can_migrant_generation", "peer_reviewed", "language",
  "year_mid_fct", "period_narrow"
)
numeric.cols <- c("year_mid")
fe_terms.main <- c("treatment_group")

housing.meta.by.ground <- split(housing.meta, as.character(housing.meta$ground_abbr))
housing.meta.scope.data <- c(list(overall = housing.meta), housing.meta.by.ground)
model.types <- c("pet", "peese")

# Helper functions ####
run_uwls_sensitivity <- function(
    scope_data,
    fe_terms,
    reg_candidates_override_fn,
    interaction_terms_fn,
    tag
) {
  housing.meta.design.list <- purrr::imap(
    scope_data,
    function(df, scope_name) {
      reg_candidates_override <- if (is.null(reg_candidates_override_fn)) {
        NULL
      } else {
        reg_candidates_override_fn(scope_name)
      }
      interaction_terms <- if (is.null(interaction_terms_fn)) {
        character(0)
      } else {
        interaction_terms_fn(scope_name)
      }

      purrr::map(
        model.types,
        function(model_type) {
          build_uwls_design(
            data = df,
            scope_name = scope_name,
            model_type = model_type,
            factor_cols = factor.cols,
            numeric_cols = numeric.cols,
            fe_terms = fe_terms,
            reg_candidates_override = reg_candidates_override,
            interaction_terms = interaction_terms
          )
        }
      )
    }
  )

  housing.meta.design <- purrr::map_dfr(
    housing.meta.design.list,
    function(scope_list) purrr::map_dfr(scope_list, "design")
  )

  housing.meta.fit.list <- purrr::imap(
    housing.meta.design.list,
    function(scope_list, scope_name) {
      purrr::map(scope_list, function(scope_obj) fit_uwls_scope(scope_obj, factor.cols))
    }
  )

  housing.meta.uwls.table.all <- purrr::map_dfr(
    housing.meta.fit.list,
    function(scope_list) purrr::map_dfr(scope_list, "table")
  )

  housing.meta.uwls.fit.status <- purrr::map_dfr(
    housing.meta.fit.list,
    function(scope_list) purrr::map_dfr(scope_list, "fit_status")
  )

  scope_names <- names(scope_data)
  housing.meta.uwls.choice <- purrr::map_dfr(
    scope_names,
    function(scope_name) {
      pet_status <- housing.meta.uwls.fit.status %>%
        filter(scope == scope_name, model_type == "pet") %>%
        pull(fit_status)
      peese_status <- housing.meta.uwls.fit.status %>%
        filter(scope == scope_name, model_type == "peese") %>%
        pull(fit_status)

      pet_p <- housing.meta.uwls.table.all %>%
        filter(model_scope == scope_name, model_type == "pet", term == "prratio_ln_se") %>%
        pull(p_value)

      pet_status <- if (length(pet_status) == 0) "missing" else pet_status[[1]]
      peese_status <- if (length(peese_status) == 0) "missing" else peese_status[[1]]
      pet_p <- if (length(pet_p) == 0) NA_real_ else pet_p[[1]]

      selected <- dplyr::case_when(
        pet_status == "ok" & is.finite(pet_p) & pet_p < 0.10 ~ "peese",
        pet_status == "ok" ~ "pet",
        pet_status != "ok" & peese_status == "ok" ~ "peese",
        TRUE ~ "none"
      )

      tibble(
        scope = scope_name,
        pet_p_value = pet_p,
        pet_status = pet_status,
        peese_status = peese_status,
        selected_model = selected
      )
    }
  )

  housing.meta.uwls.table.selected <- housing.meta.uwls.table.all %>%
    left_join(
      housing.meta.uwls.choice %>% select(scope, selected_model),
      by = c("model_scope" = "scope")
    ) %>%
    filter(model_type == selected_model)

  housing.meta.uwls.table.selected <- housing.meta.uwls.table.selected %>%
    mutate(
      model_scope = if_else(
        model_scope == "overall",
        "overall",
        paste0("ground_", model_scope)
      ),
      ground_abbr = if_else(model_scope == "overall", NA_character_, gsub("^ground_", "", model_scope))
    )

  housing.meta.uwls.table.pet <- housing.meta.uwls.table.all %>%
    filter(model_type == "pet")

  housing.meta.uwls.table.peese <- housing.meta.uwls.table.all %>%
    filter(model_type == "peese")

  list(
    tag = tag,
    design = housing.meta.design,
    fit_status = housing.meta.uwls.fit.status,
    choice = housing.meta.uwls.choice,
    tables = list(
      selected = housing.meta.uwls.table.selected,
      pet = housing.meta.uwls.table.pet,
      peese = housing.meta.uwls.table.peese
    )
  )
}

# Model runs ####
housing.meta.design.list <- purrr::imap(
  housing.meta.scope.data,
  function(df, scope_name) {
    purrr::map(
      model.types,
      function(model_type) {
        build_uwls_design(
          data = df,
          scope_name = scope_name,
          model_type = model_type,
          factor_cols = factor.cols,
          numeric_cols = numeric.cols,
          fe_terms = fe_terms.main,
          interaction_terms = character(0)
        )
      }
    )
  }
)

housing.meta.design <- purrr::map_dfr(
  housing.meta.design.list,
  function(scope_list) purrr::map_dfr(scope_list, "design")
)

housing.meta.fit.list <- purrr::imap(
  housing.meta.design.list,
  function(scope_list, scope_name) {
    purrr::map(scope_list, function(scope_obj) fit_uwls_scope(scope_obj, factor.cols))
  }
)

housing.meta.uwls.table.all <- purrr::map_dfr(
  housing.meta.fit.list,
  function(scope_list) purrr::map_dfr(scope_list, "table")
)

housing.meta.uwls.fit.status <- purrr::map_dfr(
  housing.meta.fit.list,
  function(scope_list) purrr::map_dfr(scope_list, "fit_status")
)

scope.names <- names(housing.meta.scope.data)
housing.meta.uwls.choice <- purrr::map_dfr(
  scope.names,
  function(scope_name) {
    pet_status <- housing.meta.uwls.fit.status %>%
      filter(scope == scope_name, model_type == "pet") %>%
      pull(fit_status)
    peese_status <- housing.meta.uwls.fit.status %>%
      filter(scope == scope_name, model_type == "peese") %>%
      pull(fit_status)

    pet_p <- housing.meta.uwls.table.all %>%
      filter(model_scope == scope_name, model_type == "pet", term == "prratio_ln_se") %>%
      pull(p_value)

    pet_status <- if (length(pet_status) == 0) "missing" else pet_status[[1]]
    peese_status <- if (length(peese_status) == 0) "missing" else peese_status[[1]]
    pet_p <- if (length(pet_p) == 0) NA_real_ else pet_p[[1]]

    selected <- dplyr::case_when(
      pet_status == "ok" & is.finite(pet_p) & pet_p < 0.10 ~ "peese",
      pet_status == "ok" ~ "pet",
      pet_status != "ok" & peese_status == "ok" ~ "peese",
      TRUE ~ "none"
    )

    tibble(
      scope = scope_name,
      pet_p_value = pet_p,
      pet_status = pet_status,
      peese_status = peese_status,
      selected_model = selected
    )
  }
)

housing.meta.uwls.table.selected <- housing.meta.uwls.table.all %>%
  left_join(
    housing.meta.uwls.choice %>% select(scope, selected_model),
    by = c("model_scope" = "scope")
  ) %>%
  filter(model_type == selected_model)

housing.meta.uwls.table.selected <- housing.meta.uwls.table.selected %>%
  mutate(
    model_scope = if_else(
      model_scope == "overall",
      "overall",
      paste0("ground_", model_scope)
    ),
    ground_abbr = if_else(model_scope == "overall", NA_character_, gsub("^ground_", "", model_scope))
  )

housing.meta.uwls.table.pet <- housing.meta.uwls.table.all %>%
  filter(model_type == "pet")

housing.meta.uwls.table.peese <- housing.meta.uwls.table.all %>%
  filter(model_type == "peese")

# Outputs ####
housing.meta.uwls.models <- list(
  design = housing.meta.design,
  fit_status = housing.meta.uwls.fit.status,
  pet_peese_choice = housing.meta.uwls.choice,
  tables = list(
    selected = housing.meta.uwls.table.selected,
    pet = housing.meta.uwls.table.pet,
    peese = housing.meta.uwls.table.peese
  )
)

write_uwls_output_bundle(
  uwls_obj = housing.meta.uwls.models,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$main),
  tag = uwls.specs$main$tag
)

# Sensitivity runs ####
reg_candidates_drop_year <- function(scope_name) {
  candidates <- uwls_candidates_for_scope(scope_name)
  candidates[candidates != "year_mid"]
}

reg_candidates_country <- function(scope_name) {
  candidates <- uwls_candidates_for_scope(scope_name)
  unique(candidates[candidates != "region_agg"])
}

reg_candidates_country_drop_year <- function(scope_name) {
  candidates <- reg_candidates_country(scope_name)
  candidates[candidates != "year_mid"]
}

uwls.sens.ctr <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group", "country"),
  reg_candidates_override_fn = reg_candidates_country,
  interaction_terms_fn = NULL,
  tag = uwls.specs$sens_ctr$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.ctr,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$sens_ctr),
  tag = uwls.specs$sens_ctr$tag
)

uwls.sens.ymf.region <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group", "year_mid_fct"),
  reg_candidates_override_fn = reg_candidates_drop_year,
  interaction_terms_fn = NULL,
  tag = uwls.specs$sens_ymf_region$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.ymf.region,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$sens_ymf_region),
  tag = uwls.specs$sens_ymf_region$tag
)

uwls.sens.ymf.ctr <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group", "year_mid_fct", "country"),
  reg_candidates_override_fn = reg_candidates_country_drop_year,
  interaction_terms_fn = NULL,
  tag = uwls.specs$sens_ymf_ctr$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.ymf.ctr,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$sens_ymf_ctr),
  tag = uwls.specs$sens_ymf_ctr$tag
)

uwls.sens.period.region <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group", "period_narrow"),
  reg_candidates_override_fn = reg_candidates_drop_year,
  interaction_terms_fn = NULL,
  tag = uwls.specs$sens_period_region$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.period.region,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$sens_period_region),
  tag = uwls.specs$sens_period_region$tag
)

uwls.sens.int.region.tg <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = character(0),
  reg_candidates_override_fn = NULL,
  interaction_terms_fn = function(scope_name) "region_agg:treatment_group",
  tag = uwls.specs$int_region_tg$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.int.region.tg,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$int_region_tg),
  tag = uwls.specs$int_region_tg$tag
)

uwls.sens.int.region.year <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group"),
  reg_candidates_override_fn = NULL,
  interaction_terms_fn = function(scope_name) "region_agg:year_mid",
  tag = uwls.specs$int_region_year$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.int.region.year,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$int_region_year),
  tag = uwls.specs$int_region_year$tag
)

uwls.sens.int.region.period <- run_uwls_sensitivity(
  scope_data = housing.meta.scope.data,
  fe_terms = c("treatment_group"),
  reg_candidates_override_fn = function(scope_name) {
    cands <- uwls_candidates_for_scope(scope_name)
    cands[cands == "year_mid"] <- "period_narrow"
    unique(cands)
  },
  interaction_terms_fn = function(scope_name) "region_agg:period_narrow",
  tag = uwls.specs$int_region_period$tag
)
write_uwls_output_bundle(
  uwls_obj = uwls.sens.int.region.period,
  output_dir = reg_output_dir_for_spec(dir.uwls.dirs, uwls.specs$int_region_period),
  tag = uwls.specs$int_region_period$tag
)

export_uwls_reporting_workbooks(
  dir_map = dir.uwls.dirs,
  specs = uwls.specs,
  uwls_objects = list(
    main = housing.meta.uwls.models,
    sens_ctr = uwls.sens.ctr,
    sens_ymf_region = uwls.sens.ymf.region,
    sens_ymf_ctr = uwls.sens.ymf.ctr,
    sens_period_region = uwls.sens.period.region,
    int_region_tg = uwls.sens.int.region.tg,
    int_region_year = uwls.sens.int.region.year,
    int_region_period = uwls.sens.int.region.period
  )
)

message("UWLS meta-regression complete.")
