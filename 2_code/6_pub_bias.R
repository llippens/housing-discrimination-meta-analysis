# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "dplyr", "purrr", "tibble", "readr",
  "metafor", "writexl", "puniform"
)

missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))
source(file.path(here::here(), "2_code", "functions", "core", "f_models.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.pub.bias <- file.path(dir.tables, "6_pub_bias")
ensure_dir(dir.tables)
ensure_dir(dir.pub.bias)

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/0_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "ground_abbr", "treatment_group", "study", "cluster", "prratio_ln", "prratio_ln_v",
  "callback_count_min", "application_count_min", "callback_count_maj", "application_count_maj"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Data subsets ####
ground_sets  <- split(housing.meta, housing.meta$ground_abbr)
tgroup_sets  <- split(housing.meta, housing.meta$treatment_group)
names(tgroup_sets) <- paste0("tg_", names(tgroup_sets))

housing.meta.sets <- c(
  list(overall = housing.meta),
  ground_sets,
  tgroup_sets
)

# Helper functions ####
safe_regtest <- function(rma_obj, method_name) {
  model_arg     <- if (method_name == "Egger") "lm" else "rma"
  predictor_arg <- if (method_name == "Peters") "ninv" else "sei"
  tryCatch(
    {
      rt <- metafor::regtest(rma_obj, model = model_arg, predictor = predictor_arg)

      # Asymmetry test: slope coefficient, not the intercept (corrected estimate).
      # Egger (lm): slope on sei is row 2 of summary(rt$fit)$coefficients.
      # Peters (rma): rt$zval[[1]] is already the slope on ninv.
      if (model_arg == "lm") {
        lm_coefs  <- summary(rt$fit)$coefficients
        asym_stat <- if (nrow(lm_coefs) >= 2) lm_coefs[2, "t value"]    else NA_real_
        asym_pval <- if (nrow(lm_coefs) >= 2) lm_coefs[2, "Pr(>|t|)"] else NA_real_
      } else {
        asym_stat <- rt$zval[[1]]
        asym_pval <- rt$pval[[1]]
      }

      tibble(
        method    = paste0("regtest_", tolower(method_name)),
        status    = "ok",
        statistic = asym_stat,
        p_value   = asym_pval,
        note      = NA_character_
      )
    },
    error = function(e) {
      tibble(
        method    = paste0("regtest_", tolower(method_name)),
        status    = "error",
        statistic = NA_real_,
        p_value   = NA_real_,
        note      = conditionMessage(e)
      )
    }
  )
}

diag_empty_row <- function(subset_name, method_name, status_name, note_text = NA_character_) {
  tibble(
    subset    = subset_name,
    method    = method_name,
    status    = status_name,
    statistic = NA_real_,
    p_value   = NA_real_,
    note      = note_text
  )
}

extract_named_numeric <- function(obj, candidate_names) {
  if (!is.list(obj) || is.null(names(obj))) {
    return(NA_real_)
  }

  for (candidate_name in candidate_names) {
    if (!(candidate_name %in% names(obj))) {
      next
    }

    value_obj <- obj[[candidate_name]]
    if (!is.atomic(value_obj) || !is.numeric(value_obj)) {
      next
    }

    finite_values <- value_obj[is.finite(value_obj)]
    if (length(finite_values) > 0) {
      return(finite_values[[1]])
    }
  }

  return(NA_real_)
}

safe_pcurve <- function(rma_obj, subset_name) {
  yi_vec <- as.numeric(rma_obj$yi)
  vi_vec <- as.numeric(rma_obj$vi)
  z_vec  <- yi_vec / sqrt(vi_vec)
  p_vec  <- 2 * pnorm(abs(z_vec), lower.tail = FALSE)

  expected_pos <- subset_name %in% c("gen", "tg_Female")
  sig_idx <- which(p_vec < 0.05 & if (expected_pos) yi_vec > 0 else yi_vec < 0)
  k_sig   <- length(sig_idx)

  if (k_sig < 3) {
    return(diag_empty_row(subset_name, "pcurve", "insufficient_sig_k",
      paste0("Fewer than 3 significant effects in expected direction (k_sig = ", k_sig, ")")))
  }

  pp_vals <- p_vec[sig_idx] / 0.05
  z_test  <- sum(qnorm(pp_vals)) / sqrt(k_sig)
  p_test  <- pnorm(z_test)

  diag_empty_row(subset_name, "pcurve", "ok") %>%
    mutate(
      statistic = z_test,
      p_value   = p_test,
      note      = paste0("k_sig = ", k_sig, "; side = ", if (expected_pos) "right" else "left")
    )
}

safe_puniform <- function(df, subset_name) {
  if (!requireNamespace("puniform", quietly = TRUE)) {
    return(diag_empty_row(subset_name, "puniform", "package_missing", "Package puniform not installed"))
  }

  # gen/Female: right tail; all others: left.
  side_arg <- if (subset_name %in% c("gen", "tg_Female")) "right" else "left"

  tryCatch(
    {
      pu_obj <- puniform::puni_star(yi = df$prratio_ln, vi = df$prratio_ln_v, side = side_arg, method = "ML")
      diag_empty_row(subset_name, "puniform", "ok") %>%
        mutate(
          statistic = extract_named_numeric(pu_obj, c("L.0", "z", "zval", "statistic")),
          p_value   = extract_named_numeric(pu_obj, c("pval.0", "p", "pval", "p.value")),
          note      = paste0("side = ", side_arg)
        )
    },
    error = function(e) {
      diag_empty_row(subset_name, "puniform", "error", conditionMessage(e))
    }
  )
}

safe_selmodel <- function(rma_obj, subset_name) {
  tryCatch(
    {
      sm_obj <- metafor::selmodel(rma_obj, type = "stepfun", steps = 0.025)
      diag_empty_row(subset_name, "selmodel", "ok") %>%
        mutate(
          statistic = extract_named_numeric(sm_obj, c("LRT")),
          p_value   = extract_named_numeric(sm_obj, c("LRTp"))
        )
    },
    error = function(e) {
      diag_empty_row(subset_name, "selmodel", "error", conditionMessage(e))
    }
  )
}

housing.meta.bias.results <- purrr::imap(
  housing.meta.sets,
  function(df, subset_name) {
    if (nrow(df) < 10) {
      insuff_note <- "Subset has fewer than 10 effects"
      return(list(
        regtest = bind_rows(
          diag_empty_row(subset_name, "regtest_egger", "insufficient_k", insuff_note),
          diag_empty_row(subset_name, "regtest_peters", "insufficient_k", insuff_note)
        ),
        pcurve = diag_empty_row(subset_name, "pcurve", "insufficient_k", insuff_note),
        puniform = diag_empty_row(subset_name, "puniform", "insufficient_k", insuff_note),
        selmodel = diag_empty_row(subset_name, "selmodel", "insufficient_k", insuff_note)
      ))
    }

    rma_uni <- tryCatch(
      metafor::rma(
        yi = prratio_ln, vi = prratio_ln_v,
        ni = application_count_min + application_count_maj,
        data = df, method = "REML"
      ),
      error = function(e) NULL
    )

    if (is.null(rma_uni)) {
      rma_note <- "rma.uni failed for subset"
      return(list(
        regtest = bind_rows(
          diag_empty_row(subset_name, "regtest_egger", "rma_uni_failed", rma_note),
          diag_empty_row(subset_name, "regtest_peters", "rma_uni_failed", rma_note)
        ),
        pcurve = diag_empty_row(subset_name, "pcurve", "rma_uni_failed", rma_note),
        puniform = diag_empty_row(subset_name, "puniform", "rma_uni_failed", rma_note),
        selmodel = diag_empty_row(subset_name, "selmodel", "rma_uni_failed", rma_note)
      ))
    }

    regtest_tbl <- bind_rows(
      safe_regtest(rma_uni, "Egger"),
      safe_regtest(rma_uni, "Peters")
    ) %>%
      mutate(subset = subset_name) %>%
      select(subset, method, status, statistic, p_value, note)

    list(
      regtest = regtest_tbl,
      pcurve = safe_pcurve(rma_uni, subset_name),
      puniform = safe_puniform(df, subset_name),
      selmodel = safe_selmodel(rma_uni, subset_name)
    )
  }
)

housing.meta.bias <- list(
  regtest  = dplyr::bind_rows(purrr::map(housing.meta.bias.results, "regtest")),
  pcurve   = dplyr::bind_rows(purrr::map(housing.meta.bias.results, "pcurve")),
  puniform = dplyr::bind_rows(purrr::map(housing.meta.bias.results, "puniform")),
  selmodel = dplyr::bind_rows(purrr::map(housing.meta.bias.results, "selmodel"))
)

# Outputs ####
saveRDS(housing.meta.bias, file = file.path(dir.pub.bias, "housing_meta_pub_bias.rds"))
writexl::write_xlsx(
  x = housing.meta.bias,
  path = file.path(dir.pub.bias, "meta_housing_pub_bias.xlsx")
)

message("Publication-bias diagnostics complete.")
