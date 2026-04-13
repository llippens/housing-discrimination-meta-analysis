report_safe_read_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

report_pluck_first <- function(x, candidates, required = FALSE) {
  if (is.null(x)) {
    if (required) {
      stop("Missing required object names: ", paste(candidates, collapse = ", "), call. = FALSE)
    }
    return(NULL)
  }

  for (candidate in candidates) {
    if (!is.null(x[[candidate]])) {
      return(x[[candidate]])
    }
  }

  if (required) {
    stop("Missing required object names: ", paste(candidates, collapse = ", "), call. = FALSE)
  }

  NULL
}

report_fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "")
}

report_stars <- function(p) {
  ifelse(
    !is.finite(p), "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  )
}

report_fmt_freq_cell <- function(estimate, std_error, p_value, digits = 3) {
  est_chr <- report_fmt_num(estimate, digits = digits)
  se_chr <- report_fmt_num(std_error, digits = digits)
  stars_chr <- report_stars(p_value)
  ifelse(est_chr == "" | se_chr == "", "", paste0(est_chr, stars_chr, " (", se_chr, ")"))
}

report_fmt_bayes_cell <- function(estimate, est_error, digits = 3) {
  est_chr <- report_fmt_num(estimate, digits = digits)
  ee_chr <- report_fmt_num(est_error, digits = digits)
  ifelse(est_chr == "" | ee_chr == "", "", paste0(est_chr, " (", ee_chr, ")"))
}

report_fmt_ci_cell <- function(estimate, low, high, digits = 3) {
  est_chr <- report_fmt_num(estimate, digits = digits)
  low_chr <- report_fmt_num(low, digits = digits)
  high_chr <- report_fmt_num(high, digits = digits)
  ifelse(
    est_chr == "" | low_chr == "" | high_chr == "",
    "",
    paste0(est_chr, " [", low_chr, ", ", high_chr, "]")
  )
}

report_key_moderator_prefixes <- function() {
  c(
    "intrcpt", "Intercept", "(Intercept)",
    "region_agg", "country",
    "year_mid",
    "period_narrow",
    "agent_type",
    "audit_type",
    "can_gender", "can_education", "can_employment", "can_migrant_generation",
    "prratio_ln_se", "prratio_ln_v"
  )
}

report_split_key_other <- function(tbl, exclude_vars = character(0)) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(list(key_tbl = tbl, other_vars_str = ""))
  }
  key_pfx <- report_key_moderator_prefixes()
  is_key <- function(term) {
    if (is.na(term)) return(FALSE)
    if (startsWith(term, "year_mid_fct")) return(FALSE)
    any(startsWith(term, key_pfx)) | grepl(":", term, fixed = TRUE)
  }
  key_mask <- vapply(as.character(tbl$term), is_key, logical(1))
  other_tbl <- tbl[!key_mask, ]
  other_tbl <- other_tbl[!is.na(other_tbl$term), ]
  all_other <- c("country", "year_mid_fct", "treatment_group", "matched_design", "info_binary",
                 "callback_type", "airbnb", "peer_reviewed", "language")
  all_other <- all_other[!all_other %in% exclude_vars]
  present <- vapply(
    all_other,
    function(v) any(startsWith(as.character(other_tbl$term), v), na.rm = TRUE),
    logical(1)
  )
  list(
    key_tbl        = tbl[key_mask, ],
    other_vars_str = paste(all_other[present], collapse = "; ")
  )
}

report_per_scope_other_rows <- function(tbl, scope_col = "model_scope",
                                        exclude_vars = character(0)) {
  empty_rows <- tibble::tibble(
    term = character(0), model_scope = character(0), cell = character(0)
  )
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(list(key_tbl = tbl, other_rows = empty_rows))
  }
  scopes <- unique(as.character(tbl[[scope_col]]))
  results <- lapply(scopes, function(ms) {
    tbl_ms <- tbl[as.character(tbl[[scope_col]]) == ms, ]
    sp <- report_split_key_other(tbl_ms, exclude_vars = exclude_vars)
    other_row <- if (nzchar(sp$other_vars_str)) {
      tibble::tibble(term = "Other controls", model_scope = ms, cell = sp$other_vars_str)
    } else {
      empty_rows
    }
    list(key_tbl = sp$key_tbl, other_row = other_row)
  })
  list(
    key_tbl    = dplyr::bind_rows(lapply(results, `[[`, "key_tbl")),
    other_rows = dplyr::bind_rows(lapply(results, `[[`, "other_row"))
  )
}

report_make_model_scope <- function(scope) {
  ifelse(scope == "overall", "overall", paste0("ground_", scope))
}

report_build_wide <- function(data, row_col, model_col, value_col, model_order = NULL) {
  if (is.null(data) || nrow(data) == 0) {
    return(tibble::tibble(term = character(0)))
  }

  x <- data %>%
    dplyr::select(
      row_value = dplyr::all_of(row_col),
      model_value = dplyr::all_of(model_col),
      cell_value = dplyr::all_of(value_col)
    ) %>%
    dplyr::mutate(
      row_value = as.character(row_value),
      model_value = as.character(model_value),
      cell_value = as.character(cell_value)
    ) %>%
    dplyr::filter(!is.na(row_value), row_value != "", !is.na(model_value), model_value != "") %>%
    dplyr::distinct(row_value, model_value, .keep_all = TRUE)

  if (nrow(x) == 0) {
    return(tibble::tibble(term = character(0)))
  }

  rows <- unique(x$row_value)
  models <- unique(x$model_value)
  if (!is.null(model_order) && length(model_order) > 0) {
    models <- c(model_order[model_order %in% models], models[!models %in% model_order])
  }

  out <- tibble::tibble(term = rows)
  for (m in models) {
    x_m <- x[x$model_value == m, c("row_value", "cell_value")]
    idx <- match(out$term, x_m$row_value)
    out[[m]] <- ifelse(is.na(idx), "", x_m$cell_value[idx])
  }

  out
}

report_bind_param_rows <- function(coef_wide, param_wide) {
  if (is.null(param_wide) || nrow(param_wide) == 0) {
    return(coef_wide)
  }

  cols_union <- unique(c(names(coef_wide), names(param_wide)))
  for (nm in setdiff(cols_union, names(coef_wide))) {
    coef_wide[[nm]] <- ""
  }
  for (nm in setdiff(cols_union, names(param_wide))) {
    param_wide[[nm]] <- ""
  }

  dplyr::bind_rows(
    coef_wide[, cols_union, drop = FALSE],
    param_wide[, cols_union, drop = FALSE]
  )
}

report_clean_sheet_names <- function(sheet_names) {
  cleaned <- gsub("[:\\\\/?*\\[\\]]", "_", sheet_names)
  cleaned <- gsub("\\s+", "_", cleaned)
  cleaned <- substr(cleaned, 1, 31)

  seen <- list()
  out <- character(length(cleaned))
  for (i in seq_along(cleaned)) {
    name_i <- cleaned[[i]]
    if (is.null(seen[[name_i]])) {
      seen[[name_i]] <- 1L
      out[[i]] <- name_i
    } else {
      seen[[name_i]] <- seen[[name_i]] + 1L
      suffix <- paste0("_", seen[[name_i]])
      keep <- 31 - nchar(suffix)
      out[[i]] <- paste0(substr(name_i, 1, keep), suffix)
    }
  }

  out
}

report_order_sheet_names <- function(sheet_names) {
  if (length(sheet_names) == 0) {
    return(sheet_names)
  }

  priority <- ifelse(grepl("^manuscript", sheet_names), 1L, 2L)
  sheet_names[order(priority, seq_along(sheet_names))]
}

report_write_workbook <- function(sheets, file_path) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package writexl is required for workbook export: ", file_path, call. = FALSE)
  }

  if (length(sheets) == 0) {
    message("Skipping workbook export (no sheets): ", file_path)
    return(invisible(FALSE))
  }

  ordered_names <- report_order_sheet_names(names(sheets))
  sheets <- sheets[ordered_names]
  names(sheets) <- report_clean_sheet_names(names(sheets))
  tmp_path <- paste0(file_path, ".tmp")
  on.exit(if (file.exists(tmp_path)) file.remove(tmp_path), add = TRUE)
  writexl::write_xlsx(x = sheets, path = tmp_path)
  if (file.exists(file_path)) {
    file.remove(file_path)
  }
  file.rename(tmp_path, file_path)
  invisible(TRUE)
}

report_extract_remra_tau2_from_obj <- function(obj) {
  if (is.null(obj)) {
    return(tibble::tibble(model_scope = character(0), tau2 = numeric(0)))
  }

  remra_tau2_value <- function(model_obj) {
    if (is.null(model_obj)) {
      return(NA_real_)
    }

    sigma_vals <- suppressWarnings(as.numeric(model_obj$sigma2))
    sigma_vals <- sigma_vals[is.finite(sigma_vals)]
    if (length(sigma_vals) > 0) {
      return(sum(sigma_vals))
    }

    tau_vals <- suppressWarnings(as.numeric(model_obj$tau2))
    tau_vals <- tau_vals[is.finite(tau_vals)]
    if (length(tau_vals) > 0) {
      return(tau_vals[[1]])
    }

    NA_real_
  }

  out <- list()

  overall_obj <- report_pluck_first(obj, c("housing.meta.model.overall", "model.overall", "overall"))
  if (!is.null(overall_obj) && !is.null(overall_obj$model)) {
    out[[length(out) + 1L]] <- tibble::tibble(
      model_scope = "overall",
      tau2 = remra_tau2_value(overall_obj$model)
    )
  }

  by_ground <- report_pluck_first(obj, c("housing.meta.models.by.ground", "models.by.ground", "by.ground"))
  if (!is.null(by_ground) && length(by_ground) > 0) {
    for (g in names(by_ground)) {
      model_obj <- by_ground[[g]]$model
      out[[length(out) + 1L]] <- tibble::tibble(
        model_scope = paste0("ground_", g),
        tau2 = remra_tau2_value(model_obj)
      )
    }
  }

  if (length(out) == 0) {
    return(tibble::tibble(model_scope = character(0), tau2 = numeric(0)))
  }

  dplyr::bind_rows(out)
}

report_extract_remra_tau2 <- function(rds_path) {
  obj <- report_safe_read_rds(rds_path)
  report_extract_remra_tau2_from_obj(obj)
}

report_extract_bhm_diagnostics_from_models <- function(models) {
  if (is.null(models) || length(models) == 0) {
    return(tibble::tibble(
      scope = character(0),
      max_rhat = numeric(0),
      min_ess_bulk = numeric(0),
      min_ess_tail = numeric(0),
      n_divergent = integer(0),
      n_treedepth_15plus = integer(0)
    ))
  }

  purrr::imap_dfr(
    models,
    function(model_obj, scope_name) {
      draw_sum <- tryCatch(
        {
          draws <- posterior::as_draws_array(model_obj)
          posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
        },
        error = function(e) NULL
      )

      max_rhat <- if (is.null(draw_sum) || !"rhat" %in% names(draw_sum)) NA_real_ else max(draw_sum$rhat, na.rm = TRUE)
      min_ess_bulk <- if (is.null(draw_sum) || !"ess_bulk" %in% names(draw_sum)) NA_real_ else min(draw_sum$ess_bulk, na.rm = TRUE)
      min_ess_tail <- if (is.null(draw_sum) || !"ess_tail" %in% names(draw_sum)) NA_real_ else min(draw_sum$ess_tail, na.rm = TRUE)

      nuts <- tryCatch(brms::nuts_params(model_obj), error = function(e) NULL)
      n_divergent <- if (is.null(nuts)) {
        NA_integer_
      } else {
        sum(nuts$Parameter == "divergent__" & nuts$Value > 0, na.rm = TRUE)
      }
      n_treedepth_15plus <- if (is.null(nuts)) {
        NA_integer_
      } else {
        sum(nuts$Parameter == "treedepth__" & nuts$Value >= 15, na.rm = TRUE)
      }

      tibble::tibble(
        scope = as.character(scope_name),
        max_rhat = as.numeric(max_rhat),
        min_ess_bulk = as.numeric(min_ess_bulk),
        min_ess_tail = as.numeric(min_ess_tail),
        n_divergent = as.integer(n_divergent),
        n_treedepth_15plus = as.integer(n_treedepth_15plus)
      )
    }
  )
}

report_extract_bhm_diagnostics <- function(rds_path) {
  models <- report_safe_read_rds(rds_path)
  report_extract_bhm_diagnostics_from_models(models)
}

report_build_freq_manuscript <- function(tbl, model_order, param_long = NULL) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    coef_wide <- tibble::tibble(term = character(0))
  } else {
    coef_long <- tbl %>%
      dplyr::filter(!is.na(term), term != "") %>%
      dplyr::mutate(cell = report_fmt_freq_cell(estimate, std_error, p_value)) %>%
      dplyr::select(term, model_scope, cell)

    coef_wide <- report_build_wide(
      data = coef_long,
      row_col = "term",
      model_col = "model_scope",
      value_col = "cell",
      model_order = model_order
    )
  }

  param_wide <- if (!is.null(param_long) && nrow(param_long) > 0) {
    report_build_wide(
      data = param_long,
      row_col = "term",
      model_col = "model_scope",
      value_col = "cell",
      model_order = model_order
    )
  } else {
    tibble::tibble(term = character(0))
  }

  report_bind_param_rows(coef_wide, param_wide)
}

report_build_bayes_manuscript <- function(tbl, model_order, include_ci = FALSE, param_long = NULL) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    coef_wide <- tibble::tibble(term = character(0))
  } else {
    coef_long <- tbl %>%
      dplyr::filter(!is.na(term), term != "") %>%
      dplyr::mutate(
        cell = if (include_ci) {
          report_fmt_ci_cell(estimate, q2.5, q97.5)
        } else {
          report_fmt_bayes_cell(estimate, est.error)
        }
      ) %>%
      dplyr::select(term, model_scope, cell)

    coef_wide <- report_build_wide(
      data = coef_long,
      row_col = "term",
      model_col = "model_scope",
      value_col = "cell",
      model_order = model_order
    )
  }

  param_wide <- if (!is.null(param_long) && nrow(param_long) > 0) {
    report_build_wide(
      data = param_long,
      row_col = "term",
      model_col = "model_scope",
      value_col = "cell",
      model_order = model_order
    )
  } else {
    tibble::tibble(term = character(0))
  }

  report_bind_param_rows(coef_wide, param_wide)
}
report_load_meta_inputs <- function(dir_meta_main, meta_main_obj = NULL, pet_peese_obj = NULL) {
  if (is.null(meta_main_obj)) {
    meta_main_obj <- report_safe_read_rds(file.path(dir_meta_main, "housing_meta_main.rds"))
  }
  if (is.null(pet_peese_obj)) {
    pet_peese_obj <- report_safe_read_rds(file.path(dir_meta_main, "housing_meta_main_pet_peese.rds"))
  }

  main_ground <- report_pluck_first(meta_main_obj, c("housing.meta.summary", "summary"))
  main_tg <- report_pluck_first(meta_main_obj, c("housing.meta.treatment.summary", "treatment.summary"))

  pet_selected <- report_pluck_first(meta_main_obj, c("housing.meta.pet.peese.selected", "pet.peese.selected"))
  if (is.null(pet_selected)) {
    pet_selected <- report_pluck_first(pet_peese_obj, c("selected", "housing.meta.pet.peese.selected"))
  }

  pet_choice <- report_pluck_first(meta_main_obj, c("housing.meta.pet.peese.choice", "pet.peese.choice"))
  if (is.null(pet_choice)) {
    pet_choice <- report_pluck_first(pet_peese_obj, c("choice", "housing.meta.pet.peese.choice"))
  }

  if (is.null(main_ground) || is.null(main_tg) || is.null(pet_selected) || is.null(pet_choice)) {
    message("Missing required meta-analysis input table(s) in RDS; skipping.")
    return(NULL)
  }

  list(
    main_ground = main_ground,
    main_tg = main_tg,
    pet_selected = pet_selected,
    pet_choice = pet_choice
  )
}

report_load_remra_spec <- function(out_dir, suffix, reg_obj = NULL) {
  if (is.null(reg_obj)) {
    reg_obj <- report_safe_read_rds(file.path(out_dir, paste0("housing_meta_reg", suffix, ".rds")))
  }

  overall_tbl <- report_pluck_first(reg_obj, c("housing.meta.table.overall", "table.overall"))
  design_tbl <- report_pluck_first(reg_obj, c("housing.meta.design", "design"))
  fit_tbl <- report_pluck_first(reg_obj, c("housing.meta.fit.status", "fit_status"))
  ground_tbl_list <- report_pluck_first(reg_obj, c("housing.meta.tables.by.ground", "tables.by.ground"))
  tau_tbl <- report_extract_remra_tau2_from_obj(reg_obj)

  if (nrow(tau_tbl) == 0) {
    tau_tbl <- report_extract_remra_tau2(file.path(out_dir, paste0("housing_meta_reg", suffix, ".rds")))
  }

  ground_tbls <- if (!is.null(ground_tbl_list) && length(ground_tbl_list) > 0) {
    purrr::map_dfr(ground_tbl_list, identity)
  } else {
    tibble::tibble()
  }

  list(
    overall_tbl = overall_tbl,
    design_tbl = design_tbl,
    fit_tbl = fit_tbl,
    tau_tbl = tau_tbl,
    ground_tbls = ground_tbls
  )
}

report_load_uwls_spec <- function(out_dir, suffix, uwls_obj = NULL) {
  if (is.null(uwls_obj)) {
    uwls_obj <- report_safe_read_rds(file.path(out_dir, paste0("housing_meta_uwls_models", suffix, ".rds")))
  }

  design_tbl <- report_pluck_first(uwls_obj, c("design", "housing.meta.design"))
  fit_tbl <- report_pluck_first(uwls_obj, c("fit_status", "housing.meta.fit.status"))
  choice_tbl <- report_pluck_first(uwls_obj, c("choice", "pet_peese_choice", "housing.meta.uwls.choice"))
  tables_obj <- report_pluck_first(uwls_obj, c("tables", "housing.meta.tables"))
  selected_tbl <- report_pluck_first(tables_obj, c("selected", "housing.meta.uwls.table.selected"))


  list(
    selected_tbl = selected_tbl,
    choice_tbl = choice_tbl,
    design_tbl = design_tbl,
    fit_tbl = fit_tbl
  )
}

report_load_bhm_spec <- function(out_dir, suffix, bayes_obj = NULL, design_tbl = NULL) {
  if (is.null(bayes_obj)) {
    bayes_obj <- list(
      models = report_safe_read_rds(file.path(out_dir, paste0("housing_meta_bayes_models", suffix, ".rds"))),
      summary = report_safe_read_rds(file.path(out_dir, paste0("housing_meta_bayes_summary", suffix, ".rds"))),
      fit_status = report_safe_read_rds(file.path(out_dir, paste0("housing_meta_bayes_fit_status", suffix, ".rds"))),
      status = report_safe_read_rds(file.path(out_dir, paste0("housing_meta_bayes_status", suffix, ".rds")))
    )
  }

  models_obj <- report_pluck_first(bayes_obj, c("models"))
  summary_tbl <- report_pluck_first(bayes_obj, c("summary"))
  fit_tbl <- report_pluck_first(bayes_obj, c("fit_status"))
  status_tbl <- report_pluck_first(bayes_obj, c("status"))

  if (is.null(design_tbl)) {
    design_tbl <- report_safe_read_rds(file.path(out_dir, paste0("housing_meta_bayes_design_check", suffix, ".rds")))
  }

  diag_tbl <- if (!is.null(models_obj)) {
    report_extract_bhm_diagnostics_from_models(models_obj)
  } else {
    report_extract_bhm_diagnostics(file.path(out_dir, paste0("housing_meta_bayes_models", suffix, ".rds")))
  }

  list(
    summary_tbl = summary_tbl,
    fit_tbl = fit_tbl,
    status_tbl = status_tbl,
    design_tbl = design_tbl,
    diag_tbl = diag_tbl
  )
}

export_meta_analysis_reporting_workbooks <- function(dir_meta_main, meta_main_obj = NULL, pet_peese_obj = NULL) {
  meta_inputs <- report_load_meta_inputs(
    dir_meta_main = dir_meta_main,
    meta_main_obj = meta_main_obj,
    pet_peese_obj = pet_peese_obj
  )

  main_ground <- meta_inputs$main_ground
  main_tg <- meta_inputs$main_tg
  pet_selected <- meta_inputs$pet_selected
  pet_choice <- meta_inputs$pet_choice

  if (is.null(main_ground) || is.null(main_tg) || is.null(pet_selected) || is.null(pet_choice)) {
    message("Skipping meta-analysis reporting workbooks: missing required input table(s).")
    return(invisible(NULL))
  }

  pet_ground <- pet_selected %>% dplyr::filter(cell_type == "ground")
  pet_ground_choice <- pet_choice %>% dplyr::filter(cell_type == "ground")
  pet_tg <- pet_selected %>% dplyr::filter(cell_type == "treatment_group")
  pet_tg_choice <- pet_choice %>% dplyr::filter(cell_type == "treatment_group")

  manuscript_ground <- main_ground %>%
    dplyr::mutate(
      row_id = as.character(ground),
      re = report_fmt_ci_cell(rr_mv_cr2, rr_mv_cr2_low, rr_mv_cr2_high),
      k = as.character(k_rows),
      N = as.character(k_studies)
    ) %>%
    dplyr::select(row_id, k, N, re) %>%
    dplyr::left_join(
      pet_ground %>%
        dplyr::mutate(row_id = as.character(ground), pet_peese = report_fmt_ci_cell(estimate_rr, conf_low_rr, conf_high_rr)) %>%
        dplyr::select(row_id, pet_peese, selected_model),
      by = "row_id"
    ) %>%
    dplyr::arrange(row_id)

  manuscript_tg <- main_tg %>%
    dplyr::mutate(
      row_id = paste0(as.character(ground_abbr), ": ", as.character(treatment_group)),
      re = report_fmt_ci_cell(rr_mv_cr2, rr_mv_cr2_low, rr_mv_cr2_high),
      k = as.character(k_rows),
      N = as.character(k_studies)
    ) %>%
    dplyr::select(row_id, k, N, re) %>%
    dplyr::left_join(
      pet_tg %>%
        dplyr::mutate(
          row_id = paste0(as.character(ground_abbr), ": ", as.character(treatment_group)),
          pet_peese = report_fmt_ci_cell(estimate_rr, conf_low_rr, conf_high_rr)
        ) %>%
        dplyr::select(row_id, pet_peese, selected_model),
      by = "row_id"
    ) %>%
    dplyr::arrange(row_id)

  ground_sheets <- list(
    raw_re_ground = main_ground,
    raw_pet_peese_selected_ground = pet_ground,
    raw_pet_peese_choice_ground = pet_ground_choice,
    manuscript_ground = manuscript_ground
  )

  tg_sheets <- list(
    raw_re_treatment_group = main_tg,
    raw_pet_peese_selected_treatment = pet_tg,
    raw_pet_peese_choice_treatment = pet_tg_choice,
    manuscript_treatment_group = manuscript_tg
  )

  report_write_workbook(ground_sheets, file.path(dir_meta_main, "meta_housing_ma_ground.xlsx"))
  report_write_workbook(tg_sheets, file.path(dir_meta_main, "meta_housing_ma_treatment_group.xlsx"))

  invisible(NULL)
}
export_remra_reporting_workbooks <- function(dir_map, specs, reg_objects = NULL) {
  overall_sheets <- list()
  ground_sheets <- list()

  for (spec in specs) {
    suffix <- reg_tag_suffix(spec$tag)
    spec_id <- spec$id
    out_dir <- reg_output_dir_for_spec(dir_map, spec)

    reg_obj <- if (is.null(reg_objects)) NULL else reg_objects[[spec_id]]
    remra_spec <- report_load_remra_spec(out_dir = out_dir, suffix = suffix, reg_obj = reg_obj)
    overall_tbl <- remra_spec$overall_tbl
    design_tbl <- remra_spec$design_tbl
    fit_tbl <- remra_spec$fit_tbl
    tau_tbl <- remra_spec$tau_tbl
    ground_tbls <- remra_spec$ground_tbls

    if (is.null(overall_tbl) || is.null(design_tbl)) {
      next
    }

    design_scope <- design_tbl %>%
      dplyr::mutate(
        model_scope = report_make_model_scope(as.character(scope)),
        k = n,
        N = n_studies
      ) %>%
      dplyr::select(model_scope, k, N)

    pseudo_r2_join <- if (!is.null(fit_tbl) && "pseudo_r2" %in% names(fit_tbl)) {
      fit_tbl %>%
        dplyr::mutate(model_scope = report_make_model_scope(as.character(scope))) %>%
        dplyr::select(model_scope, pseudo_r2)
    } else {
      design_scope %>% dplyr::mutate(pseudo_r2 = NA_real_) %>% dplyr::select(model_scope, pseudo_r2)
    }

    params <- design_scope %>%
      dplyr::left_join(tau_tbl, by = "model_scope") %>%
      dplyr::left_join(pseudo_r2_join, by = "model_scope") %>%
      dplyr::mutate(
        k = as.character(k),
        N = as.character(N),
        tau2 = report_fmt_num(tau2),
        pseudo_r2 = report_fmt_num(pseudo_r2)
      )

    # Split coefficient tables into key moderators vs other controls
    overall_tbl_main <- overall_tbl %>% dplyr::filter(model_scope == "overall")
    split_overall <- report_split_key_other(overall_tbl_main)
    ground_split  <- report_per_scope_other_rows(ground_tbls)

    other_row_overall <- if (nzchar(split_overall$other_vars_str)) {
      params %>% dplyr::filter(model_scope == "overall") %>%
        dplyr::transmute(term = "Other controls", model_scope, cell = split_overall$other_vars_str)
    } else {
      tibble::tibble(term = character(0), model_scope = character(0), cell = character(0))
    }

    overall_params <- dplyr::bind_rows(
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "tau^2", model_scope, cell = tau2),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "R^2 (pseudo)", model_scope, cell = pseudo_r2),
      other_row_overall
    )

    ground_params <- dplyr::bind_rows(
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "tau^2", model_scope, cell = tau2),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "R^2 (pseudo)", model_scope, cell = pseudo_r2),
      ground_split$other_rows
    )

    ground_scopes <- design_scope %>%
      dplyr::filter(model_scope != "overall") %>%
      dplyr::pull(model_scope)

    overall_sheets[[paste0("raw_overall_", spec_id)]] <- overall_tbl_main
    overall_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope == "overall")
    if (!is.null(fit_tbl)) {
      overall_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope == "overall")
    }
    overall_sheets[[paste0("manuscript_", spec_id)]] <- report_build_freq_manuscript(
      tbl = split_overall$key_tbl,
      model_order = c("overall"),
      param_long = overall_params
    )

    ground_sheets[[paste0("raw_by_ground_", spec_id)]] <- ground_tbls
    ground_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope != "overall")
    if (!is.null(fit_tbl)) {
      ground_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope != "overall")
    }
    ground_sheets[[paste0("manuscript_", spec_id)]] <- report_build_freq_manuscript(
      tbl = ground_split$key_tbl,
      model_order = unique(ground_scopes),
      param_long = ground_params
    )
  }

  report_write_workbook(overall_sheets, file.path(dir_map$main, "meta_housing_remra_overall.xlsx"))
  report_write_workbook(ground_sheets, file.path(dir_map$main, "meta_housing_remra_by_ground.xlsx"))

  invisible(NULL)
}

export_uwls_reporting_workbooks <- function(dir_map, specs, uwls_objects = NULL) {
  overall_sheets <- list()
  ground_sheets <- list()

  for (spec in specs) {
    suffix <- reg_tag_suffix(spec$tag)
    spec_id <- spec$id
    out_dir <- reg_output_dir_for_spec(dir_map, spec)

    uwls_obj <- if (is.null(uwls_objects)) NULL else uwls_objects[[spec_id]]
    uwls_spec <- report_load_uwls_spec(out_dir = out_dir, suffix = suffix, uwls_obj = uwls_obj)

    selected_tbl <- uwls_spec$selected_tbl
    choice_tbl <- uwls_spec$choice_tbl
    design_tbl <- uwls_spec$design_tbl
    fit_tbl <- uwls_spec$fit_tbl

    if (is.null(selected_tbl) || is.null(choice_tbl) || is.null(design_tbl)) {
      next
    }

    choice_scope <- choice_tbl %>%
      dplyr::mutate(model_scope = report_make_model_scope(as.character(scope))) %>%
      dplyr::select(model_scope, selected_model)

    design_selected <- design_tbl %>%
      dplyr::left_join(choice_tbl %>% dplyr::select(scope, selected_model), by = "scope") %>%
      dplyr::filter(model_type == selected_model) %>%
      dplyr::mutate(
        model_scope = report_make_model_scope(as.character(scope)),
        k = n,
        N = n_studies,
        fixed_effects = dplyr::if_else(
          is.na(fe_terms_used) | trimws(as.character(fe_terms_used)) == "",
          "(none)",
          trimws(as.character(fe_terms_used))
        )
      ) %>%
      dplyr::select(model_scope, k, N, fixed_effects)

    adj_r2_tbl <- if ("adj_r2" %in% names(selected_tbl)) {
      selected_tbl %>%
        dplyr::group_by(model_scope) %>%
        dplyr::summarise(adj_r2 = dplyr::first(adj_r2), .groups = "drop")
    } else {
      design_selected %>% dplyr::mutate(adj_r2 = NA_real_) %>% dplyr::select(model_scope, adj_r2)
    }

    params <- design_selected %>%
      dplyr::left_join(adj_r2_tbl, by = "model_scope") %>%
      dplyr::left_join(choice_scope, by = "model_scope") %>%
      dplyr::mutate(
        k = as.character(k),
        N = as.character(N),
        adj_r2 = report_fmt_num(adj_r2),
        selected_model = as.character(selected_model)
      )

    # Split coefficient tables into key moderators vs other controls
    overall_selected <- selected_tbl %>% dplyr::filter(model_scope == "overall")
    ground_selected  <- selected_tbl %>% dplyr::filter(model_scope != "overall")
    split_overall <- report_split_key_other(overall_selected)
    ground_split  <- report_per_scope_other_rows(ground_selected)

    other_row_overall <- if (nzchar(split_overall$other_vars_str)) {
      params %>% dplyr::filter(model_scope == "overall") %>%
        dplyr::transmute(term = "Other controls", model_scope, cell = split_overall$other_vars_str)
    } else {
      tibble::tibble(term = character(0), model_scope = character(0), cell = character(0))
    }

    overall_params <- dplyr::bind_rows(
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "Adj. R^2", model_scope, cell = adj_r2),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "Selected model", model_scope, cell = selected_model),
      params %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "Fixed effects", model_scope, cell = fixed_effects),
      other_row_overall
    )

    ground_params <- dplyr::bind_rows(
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "Adj. R^2", model_scope, cell = adj_r2),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "Selected model", model_scope, cell = selected_model),
      params %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "Fixed effects", model_scope, cell = fixed_effects),
      ground_split$other_rows
    )

    ground_scopes <- unique(ground_selected$model_scope)

    overall_sheets[[paste0("raw_selected_", spec_id)]] <- overall_selected
    overall_sheets[[paste0("raw_choice_", spec_id)]] <- choice_tbl %>% dplyr::filter(scope == "overall")
    overall_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope == "overall")
    if (!is.null(fit_tbl)) {
      overall_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope == "overall")
    }
    overall_sheets[[paste0("manuscript_", spec_id)]] <- report_build_freq_manuscript(
      tbl = split_overall$key_tbl,
      model_order = c("overall"),
      param_long = overall_params
    )

    ground_sheets[[paste0("raw_selected_", spec_id)]] <- ground_selected
    ground_sheets[[paste0("raw_choice_", spec_id)]] <- choice_tbl %>% dplyr::filter(scope != "overall")
    ground_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope != "overall")
    if (!is.null(fit_tbl)) {
      ground_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope != "overall")
    }
    ground_sheets[[paste0("manuscript_", spec_id)]] <- report_build_freq_manuscript(
      tbl = ground_split$key_tbl,
      model_order = ground_scopes,
      param_long = ground_params
    )
  }

  report_write_workbook(overall_sheets, file.path(dir_map$main, "meta_housing_uwlsmra_overall.xlsx"))
  report_write_workbook(ground_sheets, file.path(dir_map$main, "meta_housing_uwlsmra_by_ground.xlsx"))

  invisible(NULL)
}

export_bhm_reporting_workbooks <- function(dir_map, specs, bayes_objects = NULL, design_tables = NULL) {
  overall_sheets <- list()
  ground_sheets <- list()

  for (spec in specs) {
    suffix <- reg_tag_suffix(spec$tag)
    spec_id <- spec$id
    out_dir <- reg_output_dir_for_spec(dir_map, spec)

    bayes_obj <- if (is.null(bayes_objects)) NULL else bayes_objects[[spec_id]]
    design_tbl <- if (is.null(design_tables)) NULL else design_tables[[spec_id]]

    bhm_spec <- report_load_bhm_spec(
      out_dir = out_dir,
      suffix = suffix,
      bayes_obj = bayes_obj,
      design_tbl = design_tbl
    )

    summary_tbl <- bhm_spec$summary_tbl
    fit_tbl <- bhm_spec$fit_tbl
    status_tbl <- bhm_spec$status_tbl
    design_tbl <- bhm_spec$design_tbl
    diag_tbl <- bhm_spec$diag_tbl

    if (is.null(summary_tbl)) {
      next
    }

    summary_fixed <- summary_tbl
    if ("component" %in% names(summary_fixed)) {
      summary_fixed <- summary_fixed %>% dplyr::filter(component %in% c("fixef", "fixed"))
    }

    summary_fixed <- summary_fixed %>%
      dplyr::mutate(model_scope = report_make_model_scope(as.character(scope)))

    param_base <- summary_tbl %>%
      dplyr::mutate(model_scope = report_make_model_scope(as.character(scope))) %>%
      dplyr::group_by(model_scope) %>%
      dplyr::summarise(
        k = dplyr::first(n),
        N = dplyr::first(j_study),
        j_study = dplyr::first(j_study),
        j_treatment_group = dplyr::first(j_treatment_group),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        k = as.character(k),
        N = as.character(N),
        j_study = as.character(j_study),
        j_treatment_group = as.character(j_treatment_group)
      )

    # Extract random SD (tau) and Bayes R2 from summary_tbl
    has_component_col <- "component" %in% names(summary_tbl)

    rand_sd_rows <- if (has_component_col) {
      summary_tbl %>%
        dplyr::filter(component %in% c("random_sd", "sigma")) %>%
        dplyr::mutate(
          model_scope = report_make_model_scope(as.character(scope)),
          display_term = dplyr::case_when(
            term == "sd_study__Intercept"            ~ "tau_study",
            term == "sd_treatment_group__Intercept"  ~ "tau_treatment_group",
            term == "sigma"                          ~ "sigma",
            TRUE                                     ~ as.character(term)
          ),
          cell = report_fmt_ci_cell(estimate, q2.5, q97.5)
        ) %>%
        dplyr::select(model_scope, term = display_term, cell)
    } else {
      tibble::tibble(model_scope = character(0), term = character(0), cell = character(0))
    }

    bayes_r2_rows <- if (has_component_col) {
      summary_tbl %>%
        dplyr::filter(component == "bayes_r2") %>%
        dplyr::mutate(
          model_scope = report_make_model_scope(as.character(scope)),
          term        = "Bayes R^2",
          cell        = report_fmt_ci_cell(estimate, q2.5, q97.5)
        ) %>%
        dplyr::select(model_scope, term, cell)
    } else {
      tibble::tibble(model_scope = character(0), term = character(0), cell = character(0))
    }

    # Split coefficient tables into key moderators vs other controls
    summary_fixed_overall <- summary_fixed %>% dplyr::filter(model_scope == "overall")
    summary_fixed_ground  <- summary_fixed %>% dplyr::filter(model_scope != "overall")
    split_overall <- report_split_key_other(summary_fixed_overall)
    ground_split  <- report_per_scope_other_rows(summary_fixed_ground)

    other_row_overall <- if (nzchar(split_overall$other_vars_str)) {
      param_base %>% dplyr::filter(model_scope == "overall") %>%
        dplyr::transmute(term = "Other controls", model_scope, cell = split_overall$other_vars_str)
    } else {
      tibble::tibble(term = character(0), model_scope = character(0), cell = character(0))
    }

    overall_params <- dplyr::bind_rows(
      param_base %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      param_base %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      param_base %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "j_study", model_scope, cell = j_study),
      param_base %>% dplyr::filter(model_scope == "overall") %>% dplyr::transmute(term = "j_treatment_group", model_scope, cell = j_treatment_group),
      rand_sd_rows %>% dplyr::filter(model_scope == "overall"),
      bayes_r2_rows %>% dplyr::filter(model_scope == "overall"),
      other_row_overall
    )

    ground_params <- dplyr::bind_rows(
      param_base %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "k", model_scope, cell = k),
      param_base %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "N", model_scope, cell = N),
      param_base %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "j_study", model_scope, cell = j_study),
      param_base %>% dplyr::filter(model_scope != "overall") %>% dplyr::transmute(term = "j_treatment_group", model_scope, cell = j_treatment_group),
      rand_sd_rows %>% dplyr::filter(model_scope != "overall"),
      bayes_r2_rows %>% dplyr::filter(model_scope != "overall"),
      ground_split$other_rows
    )

    ground_scopes <- unique(summary_fixed$model_scope[summary_fixed$model_scope != "overall"])

    overall_sheets[[paste0("raw_summary_", spec_id)]] <- summary_tbl %>% dplyr::filter(scope == "overall")
    if (!is.null(fit_tbl)) {
      overall_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope == "overall")
    }
    if (!is.null(status_tbl)) {
      overall_sheets[[paste0("raw_status_", spec_id)]] <- status_tbl
    }
    if (!is.null(design_tbl)) {
      overall_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope == "overall")
    }
    overall_sheets[[paste0("raw_diagnostics_", spec_id)]] <- diag_tbl %>% dplyr::filter(scope == "overall")

    overall_sheets[[paste0("manuscript_", spec_id)]] <- report_build_bayes_manuscript(
      tbl = split_overall$key_tbl,
      model_order = c("overall"),
      include_ci = FALSE,
      param_long = overall_params
    )
    overall_sheets[[paste0("manuscript_ci_", spec_id)]] <- report_build_bayes_manuscript(
      tbl = split_overall$key_tbl,
      model_order = c("overall"),
      include_ci = TRUE,
      param_long = overall_params
    )

    ground_sheets[[paste0("raw_summary_", spec_id)]] <- summary_tbl %>% dplyr::filter(scope != "overall")
    if (!is.null(fit_tbl)) {
      ground_sheets[[paste0("raw_fit_", spec_id)]] <- fit_tbl %>% dplyr::filter(scope != "overall")
    }
    if (!is.null(status_tbl)) {
      ground_sheets[[paste0("raw_status_", spec_id)]] <- status_tbl
    }
    if (!is.null(design_tbl)) {
      ground_sheets[[paste0("raw_design_", spec_id)]] <- design_tbl %>% dplyr::filter(scope != "overall")
    }
    ground_sheets[[paste0("raw_diagnostics_", spec_id)]] <- diag_tbl %>% dplyr::filter(scope != "overall")

    ground_sheets[[paste0("manuscript_", spec_id)]] <- report_build_bayes_manuscript(
      tbl = ground_split$key_tbl,
      model_order = ground_scopes,
      include_ci = FALSE,
      param_long = ground_params
    )
    ground_sheets[[paste0("manuscript_ci_", spec_id)]] <- report_build_bayes_manuscript(
      tbl = ground_split$key_tbl,
      model_order = ground_scopes,
      include_ci = TRUE,
      param_long = ground_params
    )
  }

  report_write_workbook(overall_sheets, file.path(dir_map$main, "meta_housing_bhmra_overall.xlsx"))
  report_write_workbook(ground_sheets, file.path(dir_map$main, "meta_housing_bhmra_by_ground.xlsx"))

  invisible(NULL)
}
