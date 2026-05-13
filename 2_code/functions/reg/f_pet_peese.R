# PET-PEESE selector following Stanley & Doucouliagos (2014) and Reed (2023):
# switch to PEESE when the PET-fit bias-corrected effect at SE = 0 (covariates
# at sample means) is significantly non-null in `direction` at one-sided `alpha`.

# Expected direction lookup for a fitted cell. Priority:
#   1. treatment_group-specific override (tg_Female -> positive; else negative)
#   2. ground_abbr aggregate (gen -> positive; else negative)
#   3. bare scope (overall / unknown) -> negative
expected_direction_for_cell <- function(
    ground_abbr     = NA_character_,
    treatment_group = NA_character_
) {
  tg <- if (length(treatment_group) == 0) NA_character_ else treatment_group[[1]]
  gd <- if (length(ground_abbr) == 0)     NA_character_ else ground_abbr[[1]]

  if (!is.na(tg) && nzchar(tg)) {
    if (identical(tg, "tg_Female")) return("positive")
    return("negative")
  }
  if (!is.na(gd) && nzchar(gd) && identical(gd, "gen")) return("positive")
  "negative"
}

# Build the empty selector return tibble so column shape is stable across
# the success, degenerate, and error branches.
empty_pet_peese_choice <- function(
    direction      = "negative",
    alpha          = 0.10,
    selected_model = "none",
    selector_note  = NA_character_
) {
  tibble::tibble(
    bc_estimate        = NA_real_,
    bc_std_error       = NA_real_,
    bc_statistic       = NA_real_,
    bc_p_value         = NA_real_,
    bc_alpha           = as.numeric(alpha),
    bc_direction_used  = as.character(direction),
    selected_model     = as.character(selected_model),
    selector_note      = as.character(selector_note)
  )
}

# Canonical S&D / Reed selector. Caller owns `newdata` (the fitted frame of
# the PET model) and any column-name conventions.
select_pet_peese_canonical <- function(
    pet_fit,
    peese_fit,
    newdata,
    direction = c("negative", "positive", "two.sided"),
    alpha     = 0.10,
    se_var    = "prratio_ln_se",
    wts       = NULL
) {
  direction <- match.arg(direction)

  if (is.null(pet_fit) && is.null(peese_fit)) {
    return(empty_pet_peese_choice(
      direction      = direction,
      alpha          = alpha,
      selected_model = "none",
      selector_note  = "both_fits_missing"
    ))
  }
  if (is.null(pet_fit)) {
    return(empty_pet_peese_choice(
      direction      = direction,
      alpha          = alpha,
      selected_model = "peese",
      selector_note  = "pet_missing"
    ))
  }
  if (is.null(peese_fit)) {
    return(empty_pet_peese_choice(
      direction      = direction,
      alpha          = alpha,
      selected_model = "pet",
      selector_note  = "peese_missing"
    ))
  }

  if (missing(newdata) || is.null(newdata) || NROW(newdata) == 0L) {
    return(empty_pet_peese_choice(
      direction      = direction,
      alpha          = alpha,
      selected_model = "pet",
      selector_note  = "newdata_missing"
    ))
  }

  nd_zero_se <- newdata
  nd_zero_se[[se_var]] <- 0

  pred <- tryCatch(
    {
      args <- list(model = pet_fit, newdata = nd_zero_se)
      if (!is.null(wts)) args$wts <- wts
      do.call(marginaleffects::avg_predictions, args)
    },
    error = function(e) e
  )

  if (inherits(pred, "error") || is.null(pred) || NROW(pred) == 0L) {
    return(empty_pet_peese_choice(
      direction      = direction,
      alpha          = alpha,
      selected_model = "pet",
      selector_note  = if (inherits(pred, "error")) paste0("avg_predictions_failed: ", conditionMessage(pred)) else "avg_predictions_empty"
    ))
  }

  est <- suppressWarnings(as.numeric(pred$estimate[[1]]))
  se  <- suppressWarnings(as.numeric(pred$std.error[[1]]))
  zst <- suppressWarnings(as.numeric(pred$statistic[[1]]))
  pv  <- suppressWarnings(as.numeric(pred$p.value[[1]]))

  use_peese <- if (!is.finite(pv) || !is.finite(est)) {
    FALSE
  } else {
    switch(
      direction,
      "negative"  = (pv < alpha) && (est < 0),
      "positive"  = (pv < alpha) && (est > 0),
      "two.sided" = (pv < alpha)
    )
  }

  tibble::tibble(
    bc_estimate        = est,
    bc_std_error       = se,
    bc_statistic       = zst,
    bc_p_value         = pv,
    bc_alpha           = as.numeric(alpha),
    bc_direction_used  = as.character(direction),
    selected_model     = if (isTRUE(use_peese)) "peese" else "pet",
    selector_note      = NA_character_
  )
}
