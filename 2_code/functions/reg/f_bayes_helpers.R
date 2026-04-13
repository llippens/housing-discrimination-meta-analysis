get_named_level <- function(levels, term) {
  if (term %in% names(levels)) {
    return(as.integer(levels[[term]]))
  }
  NA_integer_
}

get_ngrp <- function(model, term) {
  val <- brms::ngrps(model)[[term]]
  if (is.null(val)) {
    return(NA_integer_)
  }
  as.integer(val)
}

prepare_bayes_data <- function(
    data,
    outcome_col = "prratio_ln",
    var_col = "prratio_ln_v",
    se_col = "se_prratio_ln",
    factor_cols = character(0)
) {
  data %>%
    dplyr::mutate(
      !!se_col := sqrt(.data[[var_col]]),
      dplyr::across(
        dplyr::all_of(factor_cols),
        ~ {
          if (is.factor(.x)) {
            .x
          } else {
            as.factor(.x)
          }
        }
      )
    ) %>%
    dplyr::filter(
      is.finite(.data[[outcome_col]]),
      is.finite(.data[[var_col]]),
      .data[[var_col]] > 0,
      is.finite(.data[[se_col]]),
      .data[[se_col]] > 0
    )
}

build_bayes_design <- function(
    data,
    scope_name,
    fixed_candidates,
    random_terms,
    outcome_col = "prratio_ln",
    var_col = "prratio_ln_v",
    se_col = "se_prratio_ln",
    factor_cols = character(0),
    interaction_terms = character(0)
) {
  data_use <- prepare_bayes_data(
    data = data,
    outcome_col = outcome_col,
    var_col = var_col,
    se_col = se_col,
    factor_cols = factor_cols
  )

  n_rows <- nrow(data_use)
  n_studies <- if ("study" %in% names(data_use)) {
    dplyr::n_distinct(data_use$study[!is.na(data_use$study) & as.character(data_use$study) != ""])
  } else {
    0L
  }
  fixed_candidates <- fixed_candidates[fixed_candidates %in% names(data_use)]
  fixed_use <- fixed_candidates[vapply(
    fixed_candidates,
    function(v) is_informative_fixed(data_use[[v]]),
    logical(1)
  )]
  interaction_candidates <- unique(interaction_terms[!is.na(interaction_terms) & interaction_terms != ""])
  interaction_use <- select_interaction_terms(
    data = data_use,
    interaction_terms = interaction_candidates,
    base_terms = fixed_use
  )

  random_levels <- vapply(
    random_terms,
    function(v) dplyr::n_distinct(data_use[[v]][!is.na(data_use[[v]])]),
    integer(1)
  )
  random_terms_used <- random_terms[random_levels > 1]
  random_terms_excluded <- setdiff(random_terms, random_terms_used)
  random_valid <- "study" %in% random_terms_used

  if (n_rows < 2) {
    df_stats <- calc_design_df(data_use, fixed_use, random_terms_used, interaction_use, random_df_method = "count_all")
    design <- tibble::tibble(
      scope = scope_name,
      n = n_rows,
      n_studies = n_studies,
      fixed_terms_candidate = paste(fixed_candidates, collapse = "; "),
      fixed_terms_used = paste(fixed_use, collapse = "; "),
      interaction_terms_candidate = paste(interaction_candidates, collapse = "; "),
      interaction_terms_used = paste(interaction_use, collapse = "; "),
      fixed_terms_dropped_df = "",
      random_terms = paste(random_terms_used, collapse = "; "),
      random_terms_excluded = paste(random_terms_excluded, collapse = "; "),
      random_levels_study = get_named_level(random_levels, "study"),
      random_levels_country = get_named_level(random_levels, "country"),
      random_levels_region_agg = get_named_level(random_levels, "region_agg"),
      random_levels_year_mid_fct = get_named_level(random_levels, "year_mid_fct"),
      random_levels_treatment_group = get_named_level(random_levels, "treatment_group"),
      fixed_df = df_stats$fixed_df,
      interaction_df = df_stats$interaction_df,
      random_df = df_stats$random_df,
      total_df = df_stats$total_df,
      df_rule = paste0(df_stats$total_df, " <= ", n_studies),
      df_ok = FALSE,
      random_levels_ok = random_valid,
      status = "insufficient_n"
    )
    return(list(
      data = data_use,
      scope = scope_name,
      fixed_terms = fixed_use,
      interaction_terms = interaction_use,
      random_terms = random_terms_used,
      status = "insufficient_n",
      message = "Fewer than 2 rows after finite/se filtering.",
      design = design
    ))
  }

  df_stats <- calc_design_df(data_use, fixed_use, random_terms_used, interaction_use, random_df_method = "count_all")
  df_ok <- df_stats$total_df <= n_studies
  model_status <- dplyr::case_when(
    !random_valid ~ "insufficient_random_levels",
    !df_ok ~ "df_exceeded",
    TRUE ~ "ready"
  )

  design <- tibble::tibble(
    scope = scope_name,
    n = n_rows,
    n_studies = n_studies,
    fixed_terms_candidate = paste(fixed_candidates, collapse = "; "),
    fixed_terms_used = paste(fixed_use, collapse = "; "),
    interaction_terms_candidate = paste(interaction_candidates, collapse = "; "),
    interaction_terms_used = paste(interaction_use, collapse = "; "),
    fixed_terms_dropped_df = "",
    random_terms = paste(random_terms_used, collapse = "; "),
    random_terms_excluded = paste(random_terms_excluded, collapse = "; "),
    random_levels_study = get_named_level(random_levels, "study"),
    random_levels_country = get_named_level(random_levels, "country"),
    random_levels_region_agg = get_named_level(random_levels, "region_agg"),
    random_levels_year_mid_fct = get_named_level(random_levels, "year_mid_fct"),
    random_levels_treatment_group = get_named_level(random_levels, "treatment_group"),
    fixed_df = df_stats$fixed_df,
    interaction_df = df_stats$interaction_df,
    random_df = df_stats$random_df,
    total_df = df_stats$total_df,
    df_rule = paste0(df_stats$total_df, " <= ", n_studies),
    df_ok = df_ok,
    random_levels_ok = random_valid,
    status = model_status
  )

  list(
    data = data_use,
    scope = scope_name,
    fixed_terms = fixed_use,
    interaction_terms = interaction_use,
    random_terms = random_terms_used,
    status = model_status,
    message = if (model_status == "ready") "Ready to fit." else model_status,
    design = design
  )
}

build_bayes_formula <- function(
    fixed_terms,
    interaction_terms = character(0),
    random_terms,
    outcome_col = "prratio_ln",
    se_col = "se_prratio_ln",
    estimate_sigma = TRUE
) {
  fixed_terms_all <- c(fixed_terms, interaction_terms)
  fixed_part <- if (length(fixed_terms_all) > 0) {
    paste(c("1", fixed_terms_all), collapse = " + ")
  } else {
    "1"
  }
  random_part <- paste0("(1 | ", random_terms, ")", collapse = " + ")
  stats::as.formula(
    paste0(outcome_col, " | se(", se_col, ", sigma = ", as.character(estimate_sigma), ") ~ ", fixed_part, " + ", random_part)
  )
}

extract_bayes_summary <- function(model, scope_name) {
  fixed_out <- tibble::as_tibble(brms::fixef(model), rownames = "term") %>%
    dplyr::mutate(scope = scope_name, component = "fixef")

  sd_out <- tibble::as_tibble(brms::posterior_summary(model, variable = "^sd_", regex = TRUE), rownames = "term") %>%
    dplyr::mutate(scope = scope_name, component = "random_sd")

  sigma_out <- tibble::as_tibble(brms::posterior_summary(model, variable = "^sigma$", regex = TRUE), rownames = "term") %>%
    dplyr::mutate(scope = scope_name, component = "sigma")

  r2_out <- tibble::as_tibble(brms::bayes_R2(model, summary = TRUE), rownames = "term") %>%
    dplyr::mutate(scope = scope_name, component = "bayes_r2")

  dplyr::bind_rows(fixed_out, sd_out, sigma_out, r2_out) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::mutate(
      n = nobs(model),
      j_study = get_ngrp(model, "study"),
      j_treatment_group = get_ngrp(model, "treatment_group")
    )
}
