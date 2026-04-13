assert_required_cols <- function(data, required_cols, context = "data") {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required columns in %s: %s",
        context,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}


