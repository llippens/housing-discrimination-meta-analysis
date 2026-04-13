compute_binary_effects <- function(
    data,
    events_treat = "callback_count_min",
    n_treat = "application_count_min",
    events_ctrl = "callback_count_maj",
    n_ctrl = "application_count_maj"
) {
  assert_required_cols(
    data = data,
    required_cols = c(events_treat, n_treat, events_ctrl, n_ctrl),
    context = "compute_binary_effects input"
  )

  out <- data %>%
    dplyr::mutate(
      callback_count_treat_raw = .data[[events_treat]],
      application_count_treat_raw = .data[[n_treat]],
      callback_count_ctrl_raw = .data[[events_ctrl]],
      application_count_ctrl_raw = .data[[n_ctrl]]
    ) %>%
    dplyr::mutate(
      noncallback_count_treat_raw = application_count_treat_raw - callback_count_treat_raw,
      noncallback_count_ctrl_raw = application_count_ctrl_raw - callback_count_ctrl_raw,
      prrate_min = callback_count_treat_raw / application_count_treat_raw,
      prrate_maj = callback_count_ctrl_raw / application_count_ctrl_raw,
      prratio = prrate_min / prrate_maj,
      prratio_ln = log(prratio),
      prratio_ln_se = sqrt(
        (1 / callback_count_treat_raw) -
          (1 / application_count_treat_raw) +
          (1 / callback_count_ctrl_raw) -
          (1 / application_count_ctrl_raw)
      ),
      prratio_ln_v = prratio_ln_se^2,
      prratio_ln_p = 1 / prratio_ln_se,
      prratio_ln_v_inv = 1 / prratio_ln_v,
      prratio_std = prratio_ln / prratio_ln_se,
      dratio = prratio - 1,
      odds_min = callback_count_treat_raw / noncallback_count_treat_raw,
      odds_maj = callback_count_ctrl_raw / noncallback_count_ctrl_raw,
      oddsratio = odds_min / odds_maj,
      oddsratio_ln = log(oddsratio),
      oddsratio_ln_se = sqrt(
        (1 / callback_count_treat_raw) +
          (1 / noncallback_count_treat_raw) +
          (1 / callback_count_ctrl_raw) +
          (1 / noncallback_count_ctrl_raw)
      )
    )

  out
}
