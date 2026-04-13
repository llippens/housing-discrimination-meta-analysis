prepare_uwls_data <- function(data, factor_cols, numeric_cols) {
  data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(factor_cols),
        ~ {
          if (is.factor(.x)) {
            .x
          } else {
            out <- stringr::str_trim(as.character(.x))
            out[out == ""] <- NA_character_
            as.factor(out)
          }
        }
      ),
      dplyr::across(dplyr::all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x))),
      prratio_ln_se = sqrt(prratio_ln_v),
      uwls_weight = 1 / prratio_ln_v
    ) %>%
    dplyr::filter(
      is.finite(prratio_ln),
      is.finite(prratio_ln_v),
      prratio_ln_v > 0,
      is.finite(prratio_ln_se),
      prratio_ln_se > 0,
      is.finite(uwls_weight),
      uwls_weight > 0,
      !is.na(study),
      as.character(study) != ""
    )
}

uwls_candidates_for_scope <- function(scope_name) {
  terms <- c(
    "region_agg", "treatment_group", "year_mid", "info_binary", "agent_type",
    "matched_design", "audit_type", "callback_type",
    "airbnb", "can_education", "can_employment",
    "peer_reviewed", "language"
  )

  if (scope_name != "gen") {
    terms <- c(terms, "can_gender")
  }
  if (scope_name == "ero") {
    terms <- c(terms, "can_migrant_generation")
  }

  unique(terms)
}

calc_uwls_df <- function(data, reg_terms, fe_terms, pet_term, interaction_terms = character(0)) {
  reg_df <- sum(vapply(reg_terms, function(v) term_df(data[[v]]), integer(1)))
  fe_df <- sum(vapply(fe_terms, function(v) term_df(data[[v]]), integer(1)))
  pet_df <- term_df(data[[pet_term]])
  interaction_df <- 0L
  if (length(interaction_terms) > 0) {
    current_terms <- c(reg_terms, pet_term)
    current_rank <- rhs_rank(data, current_terms)
    for (term in interaction_terms) {
      new_rank <- rhs_rank(data, c(current_terms, term))
      if (is.finite(current_rank) && is.finite(new_rank) && new_rank > current_rank) {
        interaction_df <- interaction_df + as.integer(new_rank - current_rank)
        current_terms <- c(current_terms, term)
        current_rank <- new_rank
      }
    }
  }
  total_df <- 1L + reg_df + fe_df + pet_df + interaction_df
  list(
    reg_df = reg_df,
    fe_df = fe_df,
    pet_df = pet_df,
    interaction_df = interaction_df,
    total_df = total_df
  )
}

build_uwls_formula <- function(reg_terms, interaction_terms, factor_cols, pet_term, fe_terms) {
  reg_terms_all <- c(reg_terms, pet_term, interaction_terms)
  reg_text <- if (length(reg_terms_all) == 0) {
    "1"
  } else {
    paste(c("1", reg_terms_all), collapse = " + ")
  }

  fe_text <- if (length(fe_terms) > 0) {
    paste(fe_terms, collapse = " + ")
  } else {
    "0"
  }

  stats::as.formula(paste0("prratio_ln ~ ", reg_text, " | ", fe_text))
}

build_uwls_design <- function(
    data,
    scope_name,
    model_type,
    factor_cols,
    numeric_cols,
    fe_terms,
    reg_candidates_override = NULL,
    interaction_terms = character(0)
) {
  pet_term <- if (model_type == "pet") "prratio_ln_se" else "prratio_ln_v"

  data_use <- prepare_uwls_data(
    data = data,
    factor_cols = factor_cols,
    numeric_cols = numeric_cols
  )

  n_rows <- nrow(data_use)
  study_levels <- dplyr::n_distinct(data_use$study[!is.na(data_use$study)])

  reg_candidates <- if (is.null(reg_candidates_override)) {
    uwls_candidates_for_scope(scope_name)
  } else {
    reg_candidates_override
  }
  reg_candidates <- reg_candidates[reg_candidates %in% names(data_use)]
  reg_use <- reg_candidates[vapply(
    reg_candidates,
    function(v) is_informative_fixed(data_use[[v]]),
    logical(1)
  )]

  fe_candidates <- fe_terms[fe_terms %in% names(data_use)]
  fe_use <- fe_candidates[vapply(
    fe_candidates,
    function(v) is_informative_fixed(data_use[[v]]),
    logical(1)
  )]
  reg_use <- setdiff(reg_use, fe_use)
  reg_use <- unique(reg_use)

  interaction_candidates <- unique(interaction_terms[!is.na(interaction_terms) & interaction_terms != ""])
  interaction_use <- select_interaction_terms(
    data = data_use,
    interaction_terms = interaction_candidates,
    base_terms = reg_use
  )

  df_stats <- calc_uwls_df(
    data = data_use,
    reg_terms = reg_use,
    fe_terms = fe_use,
    pet_term = pet_term,
    interaction_terms = interaction_use
  )
  df_ok <- df_stats$total_df <= study_levels
  model_status <- dplyr::case_when(
    n_rows < 2 ~ "insufficient_n",
    study_levels <= 1 ~ "insufficient_random_levels",
    !df_ok ~ "df_exceeded",
    TRUE ~ "ready"
  )

  design <- tibble::tibble(
    scope = scope_name,
    model_type = model_type,
    n = n_rows,
    n_studies = study_levels,
    reg_terms_candidate = paste(reg_candidates, collapse = "; "),
    reg_terms_used = paste(reg_use, collapse = "; "),
    interaction_terms_candidate = paste(interaction_candidates, collapse = "; "),
    interaction_terms_used = paste(interaction_use, collapse = "; "),
    fe_terms_used = paste(fe_use, collapse = "; "),
    pet_term = pet_term,
    reg_df = df_stats$reg_df,
    fe_df = df_stats$fe_df,
    pet_df = df_stats$pet_df,
    interaction_df = df_stats$interaction_df,
    total_df = df_stats$total_df,
    df_rule = paste0(df_stats$total_df, " <= ", study_levels),
    df_ok = df_ok,
    study_levels = study_levels,
    status = model_status
  )

  list(
    data = data_use,
    scope = scope_name,
    model_type = model_type,
    reg_terms = reg_use,
    interaction_terms = interaction_use,
    fe_terms = fe_use,
    pet_term = pet_term,
    status = model_status,
    design = design
  )
}

ssc_strict <- fixest::ssc(
  adj               = TRUE,
  t.df              = "min",
  fixef.K           = "full",
  fixef.force_exact = TRUE
)

fit_uwls_scope <- function(scope_obj, factor_cols) {
  if (scope_obj$status != "ready") {
    table_out <- tibble::tibble(
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      adj_r2 = NA_real_,
      model_scope = scope_obj$scope,
      model_type = scope_obj$model_type,
      model_formula = NA_character_,
      error = scope_obj$status
    )
    status_out <- tibble::tibble(
      scope = scope_obj$scope,
      model_type = scope_obj$model_type,
      fit_status = "skipped",
      fit_note = scope_obj$status,
      model_formula = NA_character_
    )
    return(list(table = table_out, fit_status = status_out))
  }

  formula_mod <- build_uwls_formula(
    reg_terms = scope_obj$reg_terms,
    interaction_terms = scope_obj$interaction_terms,
    factor_cols = factor_cols,
    pet_term = scope_obj$pet_term,
    fe_terms = scope_obj$fe_terms
  )
  formula_chr <- paste(format(formula_mod), collapse = " ")

  fit_obj <- tryCatch(
    fixest::feols(
      fml = formula_mod,
      data = scope_obj$data,
      weights = scope_obj$data$uwls_weight,
      vcov = fixest::vcov_cluster("study")
    ),
    error = function(e) e
  )

  if (inherits(fit_obj, "error")) {
    table_out <- tibble::tibble(
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      adj_r2 = NA_real_,
      model_scope = scope_obj$scope,
      model_type = scope_obj$model_type,
      model_formula = formula_chr,
      error = conditionMessage(fit_obj)
    )
    status_out <- tibble::tibble(
      scope = scope_obj$scope,
      model_type = scope_obj$model_type,
      fit_status = "error",
      fit_note = conditionMessage(fit_obj),
      model_formula = formula_chr
    )
    return(list(table = table_out, fit_status = status_out))
  }

  ct <- fixest::coeftable(fit_obj, vcov = ~study, ssc = ssc_strict)
  table_out <- tibble::tibble(
    term = rownames(ct),
    estimate = unname(ct[, 1]),
    std_error = unname(ct[, 2]),
    statistic = unname(ct[, 3]),
    p_value = unname(ct[, 4])
  )

  ci <- tryCatch(fixest::confint(fit_obj, vcov = ~study, ssc = ssc_strict), error = function(e) NULL)
  if (!is.null(ci)) {
    ci_tbl <- tibble::tibble(
      term = rownames(ci),
      conf_low = ci[, 1],
      conf_high = ci[, 2]
    )
    table_out <- table_out %>%
      dplyr::left_join(ci_tbl, by = "term")
  } else {
    table_out <- table_out %>%
      dplyr::mutate(conf_low = NA_real_, conf_high = NA_real_)
  }

  table_out <- table_out %>%
    dplyr::mutate(
      adj_r2 = tryCatch(as.numeric(fixest::r2(fit_obj, "ar2")), error = function(e) NA_real_),
      model_scope = scope_obj$scope,
      model_type = scope_obj$model_type,
      model_formula = formula_chr,
      error = NA_character_
    )

  status_out <- tibble::tibble(
    scope = scope_obj$scope,
    model_type = scope_obj$model_type,
    fit_status = "ok",
    fit_note = NA_character_,
    model_formula = formula_chr
  )

  list(model = fit_obj, table = table_out, fit_status = status_out)
}
