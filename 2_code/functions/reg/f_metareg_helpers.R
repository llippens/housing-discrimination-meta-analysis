infer_random_terms <- function(data) {
  if (!("cluster" %in% names(data))) {
    return("study")
  }

  study_chr <- as.character(data$study)
  cluster_chr <- as.character(data$cluster)
  valid_idx <- !is.na(study_chr) & study_chr != "" & !is.na(cluster_chr) & cluster_chr != ""

  if (!any(valid_idx)) {
    return("study")
  }

  cluster_map <- split(cluster_chr[valid_idx], study_chr[valid_idx])
  has_within_study_cluster_variation <- any(vapply(
    cluster_map,
    function(x) length(unique(x)) > 1,
    logical(1)
  ))

  if (has_within_study_cluster_variation) c("study", "cluster") else "study"
}

prepare_metareg_data <- function(data, factor_cols, numeric_cols) {
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
      dplyr::across(dplyr::all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    dplyr::filter(
      is.finite(prratio_ln),
      is.finite(prratio_ln_v),
      prratio_ln_v > 0,
      !is.na(study),
      as.character(study) != ""
    )
}

build_metareg_formula <- function(fixed_terms, interaction_terms = character(0), factor_cols) {
  terms_all <- c(fixed_terms, interaction_terms)
  if (length(terms_all) == 0) {
    return(~ 1)
  }

  stats::as.formula(paste("~", paste(terms_all, collapse = " + ")))
}

fixed_candidates_for_scope <- function(scope_name) {
  terms <- c(
    "region_agg", "year_mid", "info_binary", "agent_type",
    "matched_design", "audit_type", "treatment_group", "callback_type",
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

build_metareg_design <- function(
    data,
    scope_name,
    factor_cols,
    numeric_cols,
    fixed_candidates_override = NULL,
    interaction_terms = character(0)
) {
  data_use <- prepare_metareg_data(
    data = data,
    factor_cols = factor_cols,
    numeric_cols = numeric_cols
  )

  n_rows <- nrow(data_use)
  n_studies <- dplyr::n_distinct(
    data_use$study[!is.na(data_use$study) & as.character(data_use$study) != ""]
  )
  fixed_candidates <- if (is.null(fixed_candidates_override)) {
    fixed_candidates_for_scope(scope_name)
  } else {
    fixed_candidates_override
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

  random_terms <- infer_random_terms(data_use)
  random_levels <- vapply(
    random_terms,
    function(v) dplyr::n_distinct(data_use[[v]][!is.na(data_use[[v]])]),
    integer(1)
  )
  random_levels_study <- if ("study" %in% names(random_levels)) random_levels[["study"]] else NA_integer_
  random_levels_cluster <- if ("cluster" %in% names(random_levels)) random_levels[["cluster"]] else NA_integer_
  random_valid <- is.finite(random_levels_study) && random_levels_study > 1

  if (n_rows < 2) {
    df_stats <- calc_design_df(data_use, fixed_use, random_terms, interaction_use, random_df_method = "study_only")
    design <- tibble::tibble(
      scope = scope_name,
      n = n_rows,
      n_studies = n_studies,
      fixed_terms_candidate = paste(fixed_candidates, collapse = "; "),
      fixed_terms_used = paste(fixed_use, collapse = "; "),
      interaction_terms_candidate = paste(interaction_candidates, collapse = "; "),
      interaction_terms_used = paste(interaction_use, collapse = "; "),
      fixed_terms_dropped_df = "",
      random_terms = paste(random_terms, collapse = "; "),
      random_levels_study = random_levels_study,
      random_levels_cluster = random_levels_cluster,
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
      status = "insufficient_n",
      design = design
    ))
  }

  df_stats <- calc_design_df(data_use, fixed_use, random_terms, interaction_use, random_df_method = "study_only")
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
    random_terms = paste(random_terms, collapse = "; "),
    random_levels_study = random_levels_study,
    random_levels_cluster = random_levels_cluster,
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
    status = model_status,
    design = design
  )
}

fit_metareg_scope <- function(scope_obj, scope_name, factor_cols) {
  if (scope_obj$status != "ready") {
    table_out <- tibble::tibble(
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      model_scope = scope_name,
      model_formula = NA_character_,
      error = scope_obj$status
    )
    status_out <- tibble::tibble(
      scope = scope_name,
      fit_status = "skipped",
      fit_note = scope_obj$status,
      model_formula = NA_character_
    )
    return(list(table = table_out, fit_status = status_out))
  }

  formula_mod <- build_metareg_formula(
    fixed_terms = scope_obj$fixed_terms,
    interaction_terms = scope_obj$interaction_terms,
    factor_cols = factor_cols
  )
  formula_chr <- paste(format(formula_mod), collapse = " ")

  fit_obj <- tryCatch(
    fit_rma_mv_mod(
      data = scope_obj$data,
      mods = formula_mod
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
      df = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      model_scope = scope_name,
      model_formula = formula_chr,
      error = conditionMessage(fit_obj)
    )
    status_out <- tibble::tibble(
      scope = scope_name,
      fit_status = "error",
      fit_note = conditionMessage(fit_obj),
      model_formula = formula_chr
    )
    return(list(table = table_out, fit_status = status_out))
  }

  cr2_obj <- tryCatch(
    tidy_rma_mv_cr2(fit_obj, study = scope_obj$data$study),
    error = function(e) e
  )

  if (inherits(cr2_obj, "error")) {
    table_out <- tibble::tibble(
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      model_scope = scope_name,
      model_formula = formula_chr,
      error = conditionMessage(cr2_obj)
    )
    status_out <- tibble::tibble(
      scope = scope_name,
      fit_status = "error",
      fit_note = conditionMessage(cr2_obj),
      model_formula = formula_chr
    )
    return(list(table = table_out, fit_status = status_out))
  }

  table_out <- cr2_obj %>%
    dplyr::mutate(
      model_scope = scope_name,
      model_formula = formula_chr,
      error = NA_character_
    )

  status_out <- tibble::tibble(
    scope = scope_name,
    fit_status = "ok",
    fit_note = NA_character_,
    model_formula = formula_chr
  )

  list(
    model = fit_obj,
    table = table_out,
    fit_status = status_out
  )
}

fit_metareg_null <- function(data, outcome_col = "prratio_ln", v_col = "prratio_ln_v") {
  tryCatch({
    metafor::rma.mv(
      yi = data[[outcome_col]],
      V = data[[v_col]],
      random = infer_random_formula(data),
      method = "REML",
      test = "t",
      data = data
    )
  }, error = function(e) NULL)
}

extract_tau2_sum <- function(fit) {
  if (is.null(fit)) NA_real_ else sum(fit$sigma2, na.rm = TRUE)
}

add_pseudo_r2_to_fit_status <- function(design_list, fit_list, fit_status_tbl) {
  pseudo_rows <- purrr::imap_dfr(
    design_list,
    function(scope_obj, scope_name) {
      fit_obj <- fit_list[[scope_name]]
      if (scope_obj$status != "ready" || is.null(fit_obj) || is.null(fit_obj$model)) {
        return(tibble::tibble(
          scope = scope_name,
          tau2_null = NA_real_,
          tau2_model = NA_real_,
          pseudo_r2 = NA_real_
        ))
      }
      tau2_model <- extract_tau2_sum(fit_obj$model)
      null_fit   <- fit_metareg_null(scope_obj$data)
      tau2_null  <- extract_tau2_sum(null_fit)
      pseudo_r2  <- if (is.finite(tau2_null) && tau2_null > 0) {
        pmax(0, (tau2_null - tau2_model) / tau2_null)
      } else {
        NA_real_
      }
      tibble::tibble(
        scope      = scope_name,
        tau2_null  = tau2_null,
        tau2_model = tau2_model,
        pseudo_r2  = pseudo_r2
      )
    }
  )
  fit_status_tbl %>%
    dplyr::left_join(pseudo_rows, by = "scope")
}
