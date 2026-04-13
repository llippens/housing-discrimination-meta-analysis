reg_tag_suffix <- function(tag = NULL) {
  if (is.null(tag) || length(tag) == 0 || is.na(tag) || tag == "") {
    return("")
  }
  paste0("_", tag)
}

sanitize_scope_tag <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]+", "_", as.character(x)))
}

pluck_first <- function(x, candidates) {
  for (candidate in candidates) {
    if (!is.null(x[[candidate]])) {
      return(x[[candidate]])
    }
  }
  stop("Missing required object names: ", paste(candidates, collapse = ", "), call. = FALSE)
}

write_remra_output_bundle <- function(reg_obj, output_dir, tag = NULL) {
  suffix <- reg_tag_suffix(tag)
  saveRDS(reg_obj, file = file.path(output_dir, paste0("housing_meta_reg", suffix, ".rds")))
}

write_uwls_output_bundle <- function(uwls_obj, output_dir, tag = NULL) {
  suffix <- reg_tag_suffix(tag)
  saveRDS(uwls_obj, file = file.path(output_dir, paste0("housing_meta_uwls_models", suffix, ".rds")))
}

write_bayes_output_bundle <- function(bayes_obj, output_dir, tag = NULL) {
  suffix <- reg_tag_suffix(tag)

  models_obj <- pluck_first(bayes_obj, c("models"))
  summary_tbl <- pluck_first(bayes_obj, c("summary"))
  fit_status_tbl <- pluck_first(bayes_obj, c("fit_status"))
  status_tbl <- pluck_first(bayes_obj, c("status"))

  saveRDS(models_obj, file = file.path(output_dir, paste0("housing_meta_bayes_models", suffix, ".rds")))
  saveRDS(summary_tbl, file = file.path(output_dir, paste0("housing_meta_bayes_summary", suffix, ".rds")))
  saveRDS(fit_status_tbl, file = file.path(output_dir, paste0("housing_meta_bayes_fit_status", suffix, ".rds")))
  saveRDS(status_tbl, file = file.path(output_dir, paste0("housing_meta_bayes_status", suffix, ".rds")))
}
