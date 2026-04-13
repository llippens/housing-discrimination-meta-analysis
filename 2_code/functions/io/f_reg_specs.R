reg_output_dirs <- function(dir_stage) {
  dirs <- list(
    main = dir_stage,
    appendix = dir_stage,
    sensitivity = dir_stage,
    interactions = dir_stage
  )

  ensure_dir(dirs$main)

  dirs
}

reg_output_dir_for_spec <- function(dir_map, spec) {
  bucket <- spec$output_bucket
  if (!bucket %in% names(dir_map)) {
    stop("Unknown output bucket: ", bucket, call. = FALSE)
  }
  dir_map[[bucket]]
}

reg_specs_remra <- function() {
  list(
    main = list(id = "main", tag = NULL, output_bucket = "main"),
    ymf = list(id = "ymf", tag = "ymf", output_bucket = "sensitivity"),
    ctr = list(id = "ctr", tag = "ctr", output_bucket = "sensitivity"),
    int_region_tg = list(id = "int_region_tg", tag = "int_region_tg", output_bucket = "interactions"),
    int_region_year = list(id = "int_region_year", tag = "int_region_year", output_bucket = "interactions"),
    period = list(id = "period", tag = "period", output_bucket = "sensitivity"),
    int_region_period = list(id = "int_region_period", tag = "int_region_period", output_bucket = "interactions")
  )
}

reg_specs_uwls <- function() {
  list(
    main = list(id = "main", tag = NULL, output_bucket = "main"),
    sens_ctr = list(id = "sens_ctr", tag = "sens_ctr", output_bucket = "sensitivity"),
    sens_ymf_region = list(id = "sens_ymf_region", tag = "sens_ymf_region", output_bucket = "sensitivity"),
    sens_ymf_ctr = list(id = "sens_ymf_ctr", tag = "sens_ymf_ctr", output_bucket = "sensitivity"),
    int_region_tg = list(id = "int_region_tg", tag = "int_region_tg", output_bucket = "interactions"),
    int_region_year = list(id = "int_region_year", tag = "int_region_year", output_bucket = "interactions"),
    sens_period_region = list(id = "sens_period_region", tag = "sens_period_region", output_bucket = "sensitivity"),
    int_region_period = list(id = "int_region_period", tag = "int_region_period", output_bucket = "interactions")
  )
}

reg_specs_bhmra <- function() {
  list(
    main = list(id = "main", tag = NULL, output_bucket = "main"),
    ctr = list(id = "ctr", tag = "ctr", output_bucket = "sensitivity"),
    int_region_tg = list(id = "int_region_tg", tag = "int_region_tg", output_bucket = "interactions"),
    int_region_year = list(id = "int_region_year", tag = "int_region_year", output_bucket = "interactions"),
    period = list(id = "period", tag = "period", output_bucket = "sensitivity"),
    int_region_period = list(id = "int_region_period", tag = "int_region_period", output_bucket = "interactions")
  )
}
