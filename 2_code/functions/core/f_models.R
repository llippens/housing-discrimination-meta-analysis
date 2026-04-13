infer_random_formula <- function(data) {
  if (!("cluster" %in% names(data))) {
    return(~ 1 | study)
  }

  study_chr <- as.character(data$study)
  cluster_chr <- as.character(data$cluster)
  valid_idx <- !is.na(study_chr) & study_chr != "" & !is.na(cluster_chr) & cluster_chr != ""

  if (!any(valid_idx)) {
    return(~ 1 | study)
  }

  cluster_map <- split(cluster_chr[valid_idx], study_chr[valid_idx])
  has_within_study_cluster_variation <- any(vapply(
    cluster_map,
    function(x) length(unique(x)) > 1,
    logical(1)
  ))

  if (has_within_study_cluster_variation) ~ 1 | study/cluster else ~ 1 | study
}

fit_rma_mv_intercept <- function(data, yi = "prratio_ln", vi = "prratio_ln_v") {
  stopifnot(requireNamespace("metafor", quietly = TRUE))
  assert_required_cols(data, c(yi, vi, "study"))

  random_formula <- infer_random_formula(data)

  metafor::rma.mv(
    yi = data[[yi]],
    V = data[[vi]],
    random = random_formula,
    method = "REML",
    test = "t",
    data = data
  )
}

fit_rma_mv_mod <- function(data, mods, yi = "prratio_ln", vi = "prratio_ln_v") {
  stopifnot(requireNamespace("metafor", quietly = TRUE))
  assert_required_cols(data, c(yi, vi, "study"))

  random_formula <- infer_random_formula(data)

  metafor::rma.mv(
    yi = data[[yi]],
    V = data[[vi]],
    mods = mods,
    random = random_formula,
    method = "REML",
    test = "t",
    data = data
  )
}

tidy_rma_mv_cr2 <- function(model, study) {
  stopifnot(requireNamespace("clubSandwich", quietly = TRUE))

  out <- clubSandwich::coef_test(
    obj = model,
    vcov = "CR2",
    cluster = study
  )

  tibble::tibble(
    term = out$Coef,
    estimate = out$beta,
    std_error = out$SE,
    statistic = out$tstat,
    p_value = out$p_Satt,
    df = out$df_Satt,
    conf_low = estimate - stats::qt(0.975, df = df) * std_error,
    conf_high = estimate + stats::qt(0.975, df = df) * std_error
  )
}
