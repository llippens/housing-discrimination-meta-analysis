# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "dplyr", "purrr", "tibble", "stringr",
  "metafor", "clubSandwich", "fixest", "writexl"
)

missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_models.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reg_output.R"))
source(file.path(here::here(), "2_code", "functions", "io", "f_reporting_workbooks.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.meta.main <- file.path(dir.tables, "2_meta_analysis")
dir.meta.main.appendix <- file.path(dir.meta.main, "appendix")

ensure_dir(dir.tables)
ensure_dir(dir.meta.main)
ensure_dir(dir.meta.main.appendix)

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/0_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "study", "cluster", "ground", "ground_abbr", "treatment_group",
  "callback_count_min", "application_count_min", "callback_count_maj", "application_count_maj",
  "prratio_ln", "prratio_ln_v"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Helper functions ####
extract_rma_tau2_i2 <- function(rma_obj, vi_vec) {
  if (is.null(rma_obj)) return(list(tau2 = NA_real_, i2 = NA_real_))
  tau2 <- sum(rma_obj$sigma2)
  vi <- vi_vec[is.finite(vi_vec) & vi_vec > 0]
  if (length(vi) == 0) return(list(tau2 = tau2, i2 = NA_real_))
  wi <- 1 / vi
  k  <- length(wi)
  typical_v <- (k - 1) * sum(wi) / (sum(wi)^2 - sum(wi^2))
  i2 <- 100 * tau2 / (tau2 + typical_v)
  list(tau2 = tau2, i2 = max(0, min(100, i2)))
}

extract_rma_pred_interval <- function(rma_obj) {
  if (is.null(rma_obj)) return(list(pred_low = NA_real_, pred_high = NA_real_))
  pred <- tryCatch(predict(rma_obj, level = 0.95), error = function(e) NULL)
  if (is.null(pred)) return(list(pred_low = NA_real_, pred_high = NA_real_))
  list(
    pred_low  = if (is.finite(pred$cr.lb)) exp(pred$cr.lb) else NA_real_,
    pred_high = if (is.finite(pred$cr.ub)) exp(pred$cr.ub) else NA_real_
  )
}

summarise_ground <- function(df, rma_obj, rma_cr2) {
  het <- extract_rma_tau2_i2(rma_obj, df$prratio_ln_v)
  tibble::tibble(
    ground = as.character(df$ground[[1]]),
    ground_abbr = as.character(df$ground_abbr[[1]]),
    k_rows = nrow(df),
    k_studies = dplyr::n_distinct(df$study),
    rr_mv_cr2 = if (!is.null(rma_cr2)) exp(rma_cr2$estimate[rma_cr2$term == "intrcpt"][1]) else NA_real_,
    rr_mv_cr2_low = if (!is.null(rma_cr2)) exp(rma_cr2$conf_low[rma_cr2$term == "intrcpt"][1]) else NA_real_,
    rr_mv_cr2_high = if (!is.null(rma_cr2)) exp(rma_cr2$conf_high[rma_cr2$term == "intrcpt"][1]) else NA_real_,
    rr_mv_cr2_p = if (!is.null(rma_cr2)) rma_cr2$p_value[rma_cr2$term == "intrcpt"][1] else NA_real_,
    tau2 = het$tau2,
    i2 = het$i2
  )
}

prepare_pet_peese_data <- function(df) {
  df %>%
    mutate(
      study = as.character(study),
      prratio_ln = suppressWarnings(as.numeric(prratio_ln)),
      prratio_ln_v = suppressWarnings(as.numeric(prratio_ln_v)),
      prratio_ln_se = sqrt(prratio_ln_v),
      uwls_weight = 1 / prratio_ln_v
    ) %>%
    filter(
      is.finite(prratio_ln),
      is.finite(prratio_ln_v),
      prratio_ln_v > 0,
      is.finite(prratio_ln_se),
      prratio_ln_se > 0,
      is.finite(uwls_weight),
      uwls_weight > 0,
      !is.na(study),
      study != ""
    )
}

extract_coef_value <- function(mat, term, col) {
  if (is.null(mat) || !(term %in% rownames(mat))) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(mat[term, col]))
}

fit_pet_peese_cell <- function(df, model_type = c("pet", "peese")) {
  model_type <- match.arg(model_type)
  pet_term <- if (model_type == "pet") "prratio_ln_se" else "prratio_ln_v"
  formula_mod <- stats::as.formula(paste0("prratio_ln ~ 1 + ", pet_term, " | 0"))
  formula_chr <- paste(format(formula_mod), collapse = " ")

  fit_obj <- tryCatch(
    fixest::feols(
      fml = formula_mod,
      data = df,
      weights = df$uwls_weight,
      vcov = fixest::vcov_cluster("study")
    ),
    error = function(e) e
  )

  if (inherits(fit_obj, "error")) {
    return(
      tibble(
        model_type = model_type,
        fit_status = "error",
        fit_note = conditionMessage(fit_obj),
        pet_term = pet_term,
        model_formula = formula_chr,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        estimate_rr = NA_real_,
        conf_low_rr = NA_real_,
        conf_high_rr = NA_real_,
        pet_term_estimate = NA_real_,
        pet_term_p_value = NA_real_
      )
    )
  }

  ct <- tryCatch(fixest::coeftable(fit_obj), error = function(e) NULL)
  ci <- tryCatch(fixest::confint(fit_obj), error = function(e) NULL)

  estimate <- extract_coef_value(ct, "(Intercept)", 1)
  std_error <- extract_coef_value(ct, "(Intercept)", 2)
  statistic <- extract_coef_value(ct, "(Intercept)", 3)
  p_value <- extract_coef_value(ct, "(Intercept)", 4)
  conf_low <- extract_coef_value(ci, "(Intercept)", 1)
  conf_high <- extract_coef_value(ci, "(Intercept)", 2)
  if (!is.finite(conf_low) || !is.finite(conf_high)) {
    z_crit <- stats::qnorm(0.975)
    if (is.finite(estimate) && is.finite(std_error)) {
      conf_low <- estimate - z_crit * std_error
      conf_high <- estimate + z_crit * std_error
    }
  }

  tibble(
    model_type = model_type,
    fit_status = "ok",
    fit_note = NA_character_,
    pet_term = pet_term,
    model_formula = formula_chr,
    estimate = estimate,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,
    conf_low = conf_low,
    conf_high = conf_high,
    estimate_rr = if_else(is.finite(estimate), exp(estimate), NA_real_),
    conf_low_rr = if_else(is.finite(conf_low), exp(conf_low), NA_real_),
    conf_high_rr = if_else(is.finite(conf_high), exp(conf_high), NA_real_),
    pet_term_estimate = extract_coef_value(ct, pet_term, 1),
    pet_term_p_value = extract_coef_value(ct, pet_term, 4)
  )
}

choose_pet_peese <- function(model_tbl) {
  pet_status <- model_tbl %>%
    filter(model_type == "pet") %>%
    pull(fit_status)
  peese_status <- model_tbl %>%
    filter(model_type == "peese") %>%
    pull(fit_status)
  pet_p <- model_tbl %>%
    filter(model_type == "pet") %>%
    pull(pet_term_p_value)

  pet_status <- if (length(pet_status) == 0) "missing" else pet_status[[1]]
  peese_status <- if (length(peese_status) == 0) "missing" else peese_status[[1]]
  pet_p <- if (length(pet_p) == 0) NA_real_ else pet_p[[1]]

  selected_model <- dplyr::case_when(
    pet_status == "ok" & is.finite(pet_p) & pet_p < 0.10 ~ "peese",
    pet_status == "ok" ~ "pet",
    pet_status != "ok" & peese_status == "ok" ~ "peese",
    TRUE ~ "none"
  )

  tibble(
    pet_p_value = pet_p,
    pet_status = pet_status,
    peese_status = peese_status,
    selected_model = selected_model
  )
}

# Models ####
housing.meta.by.ground <- split(housing.meta, housing.meta$ground_abbr)

housing.meta.models <- purrr::imap(
  housing.meta.by.ground,
  function(df, g) {
    rma <- tryCatch(fit_rma_mv_intercept(df), error = function(e) NULL)
    cr2 <- if (!is.null(rma)) {
      tryCatch(tidy_rma_mv_cr2(rma, study = df$study), error = function(e) NULL)
    } else {
      NULL
    }
    list(data = df, rma = rma, cr2 = cr2)
  }
)

housing.meta.summary <- purrr::map_dfr(
  housing.meta.models,
  function(obj) summarise_ground(obj$data, obj$rma, obj$cr2)
) %>%
  arrange(desc(k_rows))

housing.meta.treatment.cells <- housing.meta %>%
  group_by(ground, ground_abbr, treatment_group) %>%
  summarise(k_rows = n(), k_studies = n_distinct(study), .groups = "drop") %>%
  filter(k_rows > 10, k_studies >= 2)

housing.meta.treatment.full <- purrr::pmap_dfr(
  housing.meta.treatment.cells,
  function(ground, ground_abbr, treatment_group, k_rows, k_studies) {
    df <- housing.meta %>%
      filter(
        ground_abbr == !!ground_abbr,
        treatment_group == !!treatment_group
      )
    rma <- tryCatch(fit_rma_mv_intercept(df), error = function(e) NULL)
    cr2 <- if (!is.null(rma)) {
      tryCatch(tidy_rma_mv_cr2(rma, study = df$study), error = function(e) NULL)
    } else {
      NULL
    }
    het <- extract_rma_tau2_i2(rma, df$prratio_ln_v)
    pi  <- extract_rma_pred_interval(rma)
    tibble::tibble(
      ground          = ground,
      ground_abbr     = ground_abbr,
      treatment_group = treatment_group,
      k_rows          = as.integer(k_rows),
      k_studies       = as.integer(k_studies),
      rr_mv_cr2       = if (!is.null(cr2)) exp(cr2$estimate[cr2$term == "intrcpt"][1]) else NA_real_,
      rr_mv_cr2_low   = if (!is.null(cr2)) exp(cr2$conf_low[cr2$term  == "intrcpt"][1]) else NA_real_,
      rr_mv_cr2_high  = if (!is.null(cr2)) exp(cr2$conf_high[cr2$term == "intrcpt"][1]) else NA_real_,
      rr_mv_cr2_p     = if (!is.null(cr2)) cr2$p_value[cr2$term       == "intrcpt"][1]  else NA_real_,
      tau2            = het$tau2,
      i2              = het$i2,
      pred_low        = pi$pred_low,
      pred_high       = pi$pred_high
    )
  }
)

housing.meta.treatment.summary <- housing.meta.treatment.full %>%
  dplyr::select(
    ground, ground_abbr, treatment_group, k_rows, k_studies,
    rr_mv_cr2, rr_mv_cr2_low, rr_mv_cr2_high, rr_mv_cr2_p,
    tau2, i2
  ) %>%
  arrange(ground_abbr, rr_mv_cr2)

housing.meta.pred.ground <- purrr::imap_dfr(
  housing.meta.models,
  function(obj, g) {
    pi <- extract_rma_pred_interval(obj$rma)
    tibble::tibble(
      ground      = as.character(obj$data$ground[[1]]),
      ground_abbr = as.character(g),
      k_rows      = nrow(obj$data),
      k_studies   = dplyr::n_distinct(obj$data$study),
      pred_low    = pi$pred_low,
      pred_high   = pi$pred_high
    )
  }
)

housing.meta.pred.tgroup <- housing.meta.treatment.full %>%
  dplyr::select(ground, ground_abbr, treatment_group, k_rows, k_studies, pred_low, pred_high)

housing.meta.pet.peese.cells.ground <- purrr::imap_dfr(
  housing.meta.by.ground,
  function(df, g) {
    tibble(
      cell_type = "ground",
      cell_id = paste0("ground_", g),
      ground = as.character(df$ground[[1]]),
      ground_abbr = as.character(g),
      treatment_group = NA_character_,
      k_rows = nrow(df),
      k_studies = n_distinct(df$study),
      data = list(df)
    )
  }
)

housing.meta.pet.peese.cells.tgroup <- purrr::pmap_dfr(
  housing.meta.treatment.cells,
  function(ground, ground_abbr, treatment_group, k_rows, k_studies) {
    cell_tag <- tolower(gsub("[^a-zA-Z0-9]+", "_", treatment_group))
    tibble(
      cell_type = "treatment_group",
      cell_id = paste0("tgroup_", ground_abbr, "_", cell_tag),
      ground = as.character(ground),
      ground_abbr = as.character(ground_abbr),
      treatment_group = as.character(treatment_group),
      k_rows = as.integer(k_rows),
      k_studies = as.integer(k_studies),
      data = list(
        housing.meta %>%
          filter(
            ground_abbr == !!ground_abbr,
            treatment_group == !!treatment_group
          )
      )
    )
  }
)

housing.meta.pet.peese.cells <- bind_rows(
  housing.meta.pet.peese.cells.ground,
  housing.meta.pet.peese.cells.tgroup
)

housing.meta.pet.peese.table <- purrr::pmap_dfr(
  housing.meta.pet.peese.cells,
  function(cell_type, cell_id, ground, ground_abbr, treatment_group, k_rows, k_studies, data) {
    df_use <- prepare_pet_peese_data(data)
    study_levels <- n_distinct(df_use$study)

    if (nrow(df_use) < 2 || study_levels <= 1) {
      model_tbl <- tibble(
        model_type = c("pet", "peese"),
        fit_status = "skipped",
        fit_note = if_else(nrow(df_use) < 2, "insufficient_n", "insufficient_random_levels"),
        pet_term = c("prratio_ln_se", "prratio_ln_v"),
        model_formula = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        estimate_rr = NA_real_,
        conf_low_rr = NA_real_,
        conf_high_rr = NA_real_,
        pet_term_estimate = NA_real_,
        pet_term_p_value = NA_real_
      )
    } else {
      model_tbl <- bind_rows(
        fit_pet_peese_cell(df_use, "pet"),
        fit_pet_peese_cell(df_use, "peese")
      )
    }

    choice <- choose_pet_peese(model_tbl)

    model_tbl %>%
      mutate(
        cell_type = cell_type,
        cell_id = cell_id,
        ground = ground,
        ground_abbr = ground_abbr,
        treatment_group = treatment_group,
        k_rows = as.integer(k_rows),
        k_studies = as.integer(k_studies),
        pet_p_value = choice$pet_p_value[[1]],
        pet_status = choice$pet_status[[1]],
        peese_status = choice$peese_status[[1]],
        selected_model = choice$selected_model[[1]]
      ) %>%
      select(
        cell_type, cell_id, ground, ground_abbr, treatment_group, k_rows, k_studies,
        model_type, selected_model, pet_p_value, pet_status, peese_status,
        fit_status, fit_note, pet_term, model_formula,
        estimate, estimate_rr, std_error, statistic, p_value,
        conf_low, conf_high, conf_low_rr, conf_high_rr,
        pet_term_estimate, pet_term_p_value
      )
  }
)

housing.meta.pet.peese.choice <- housing.meta.pet.peese.table %>%
  distinct(
    cell_type, cell_id, ground, ground_abbr, treatment_group, k_rows, k_studies,
    pet_p_value, pet_status, peese_status, selected_model
  ) %>%
  arrange(cell_type, ground_abbr, treatment_group)

housing.meta.pet.peese.selected <- housing.meta.pet.peese.table %>%
  filter(model_type == selected_model) %>%
  arrange(cell_type, ground_abbr, treatment_group)

# Outputs ####
housing.meta.main <- list(
  housing.meta.models              = housing.meta.models,
  housing.meta.summary             = housing.meta.summary,
  housing.meta.treatment.summary   = housing.meta.treatment.summary,
  housing.meta.pred.ground         = housing.meta.pred.ground,
  housing.meta.pred.tgroup         = housing.meta.pred.tgroup,
  housing.meta.pet.peese.table     = housing.meta.pet.peese.table,
  housing.meta.pet.peese.choice    = housing.meta.pet.peese.choice,
  housing.meta.pet.peese.selected  = housing.meta.pet.peese.selected
)

saveRDS(housing.meta.main, file = file.path(dir.meta.main, "housing_meta_main.rds"))
writexl::write_xlsx(
  list(
    summary_ground     = housing.meta.summary,
    summary_treatment  = housing.meta.treatment.summary,
    pred_ground        = housing.meta.pred.ground,
    pred_treatment     = housing.meta.pred.tgroup,
    pet_peese_table    = housing.meta.pet.peese.table,
    pet_peese_choice   = housing.meta.pet.peese.choice,
    pet_peese_selected = housing.meta.pet.peese.selected
  ),
  path = file.path(dir.meta.main, "meta_housing_meta_analysis.xlsx")
)

export_meta_analysis_reporting_workbooks(
  dir_meta_main = dir.meta.main,
  meta_main_obj = housing.meta.main
)

message("Main meta-analysis complete: ", nrow(housing.meta.summary), " ground models and ", nrow(housing.meta.treatment.summary), " treatment-group cells.")
