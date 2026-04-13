# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c("here", "dplyr", "tibble", "writexl", "modelsummary", "stringr", "purrr")
missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))

# Reproducibility ####
set.seed(8888)

# Directories ####
dir.root <- here::here()
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.descriptives <- file.path(dir.tables, "1_descriptives")
dir.descriptives.appendix <- file.path(dir.descriptives, "appendix")

ensure_dir(dir.tables)
ensure_dir(dir.descriptives)
ensure_dir(dir.descriptives.appendix)

# Inputs ####
housing.meta.path <- file.path(dir.public, "housing_meta.rds")
if (!file.exists(housing.meta.path)) {
  stop("Run 2_code/0_wrangling.R first. Missing: ", housing.meta.path, call. = FALSE)
}

housing.meta <- readRDS(housing.meta.path)

required.cols <- c(
  "application_count_min", "application_count_maj",
  "callback_count_min", "callback_count_maj",
  "prrate_min", "prrate_maj", "prratio_ln",
  "year_research_start", "year_research_end",
  "ground", "region_agg", "info_binary", "agent_type",
  "matched_design", "audit_type", "treatment_group", "callback_type", "airbnb",
  "can_gender", "can_education", "can_employment", "can_migrant_generation",
  "peer_reviewed", "language",
  "country", "country_iso3", "region",
  "period_narrow", "study"
)
assert_required_cols(housing.meta, required.cols, "housing.meta data")

# Helper functions ####
format_n <- function(x) {
  ifelse(is.na(x), "", as.character(as.integer(round(x, 0))))
}

format_num2 <- function(x) {
  ifelse(is.na(x), "", sprintf("%.2f", round(x, 2)))
}

blank_row <- function() {
  tibble(
    variable = "",
    n = "",
    mean = "",
    sd = "",
    min = "",
    median = "",
    max = "",
    percent = ""
  )
}

summarise_numeric <- function(
  df,
  var_name,
  var_label,
  require_unit_interval = FALSE,
  use_decimal_minmax = FALSE
) {
  raw_values <- suppressWarnings(as.numeric(df[[var_name]]))
  raw_values[!is.finite(raw_values)] <- NA_real_

  if (require_unit_interval) {
    valid_values <- raw_values[!is.na(raw_values)]
    if (length(valid_values) == 0 || any(valid_values < 0 | valid_values > 1)) {
      message("Skipping ", var_name, " because values are not valid rates in [0, 1].")
      return(NULL)
    }
  }

  n_valid <- sum(!is.na(raw_values))
  if (n_valid == 0) {
    return(NULL)
  }

  sd_value <- if (n_valid >= 2) stats::sd(raw_values, na.rm = TRUE) else NA_real_
  min_value <- min(raw_values, na.rm = TRUE)
  median_value <- stats::quantile(raw_values, probs = 0.50, na.rm = TRUE, names = FALSE)
  max_value <- max(raw_values, na.rm = TRUE)

  if (require_unit_interval || use_decimal_minmax) {
    min_out <- format_num2(min_value)
    median_out <- format_num2(median_value)
    max_out <- format_num2(max_value)
  } else {
    min_out <- format_n(min_value)
    median_out <- format_n(median_value)
    max_out <- format_n(max_value)
  }

  tibble(
    variable = var_label,
    n = format_n(n_valid),
    mean = format_num2(mean(raw_values, na.rm = TRUE)),
    sd = format_num2(sd_value),
    min = min_out,
    median = median_out,
    max = max_out,
    percent = ""
  )
}

# Table 1 ####
continuous.map <- tibble::tribble(
  ~var_name, ~var_label, ~require_unit_interval, ~use_decimal_minmax,
  "application_count_min", "Applications: Minority Group", FALSE, FALSE,
  "application_count_maj", "Applications: Majority Group", FALSE, FALSE,
  "callback_count_min", "Callbacks: Minority Group", FALSE, FALSE,
  "callback_count_maj", "Callbacks: Majority Group", FALSE, FALSE,
  "prrate_min", "Positive Response Rate: Minority Group", TRUE, FALSE,
  "prrate_maj", "Positive Response Rate: Majority Group", TRUE, FALSE,
  "prratio_ln", "Log Positive Response Ratio", FALSE, TRUE,
  "year_research_start", "Start Year of Data Collection", FALSE, FALSE,
  "year_research_end", "End Year of Data Collection", FALSE, FALSE
)

housing.meta.panel.a <- purrr::pmap_dfr(
  continuous.map,
  function(var_name, var_label, require_unit_interval, use_decimal_minmax) {
    summarise_numeric(
      housing.meta,
      var_name,
      var_label,
      require_unit_interval,
      use_decimal_minmax
    )
  }
)

summarise_categorical <- function(df, var_name, var_label) {
  category_counts <- df %>%
    mutate(category_value = as.character(.data[[var_name]])) %>%
    filter(!is.na(category_value), stringr::str_trim(category_value) != "") %>%
    count(category_value, name = "n") %>%
    arrange(desc(n), category_value)

  header_row <- tibble(
    variable = var_label,
    n = "",
    mean = "",
    sd = "",
    min = "",
    median = "",
    max = "",
    percent = ""
  )

  if (nrow(category_counts) == 0) {
    return(header_row)
  }

  category_rows <- category_counts %>%
    mutate(
      variable = category_value,
      percent = format_num2(100 * n / sum(n)),
      n = format_n(n),
      mean = "",
      sd = "",
      min = "",
      median = "",
      max = ""
    ) %>%
    select(variable, n, mean, sd, min, median, max, percent)

  bind_rows(header_row, category_rows)
}

summarise_ground_treatment <- function(df) {
  header_row <- tibble(
    variable = "Ground",
    n = "",
    mean = "",
    sd = "",
    min = "",
    median = "",
    max = "",
    percent = ""
  )

  ground_counts <- df %>%
    mutate(ground_value = as.character(ground)) %>%
    filter(!is.na(ground_value), stringr::str_trim(ground_value) != "") %>%
    count(ground_value, name = "n_ground") %>%
    arrange(desc(n_ground), ground_value)

  if (nrow(ground_counts) == 0) {
    return(header_row)
  }

  total_n <- sum(ground_counts$n_ground)

  ground_rows <- purrr::map_dfr(
    seq_len(nrow(ground_counts)),
    function(i) {
      ground_key <- ground_counts$ground_value[[i]]
      ground_n <- ground_counts$n_ground[[i]]

      ground_row <- tibble(
        variable = ground_key,
        n = format_n(ground_n),
        mean = "",
        sd = "",
        min = "",
        median = "",
        max = "",
        percent = format_num2(100 * ground_n / total_n)
      )

      treatment_counts <- df %>%
        mutate(
          ground_value = as.character(ground),
          treatment_value = as.character(treatment_group)
        ) %>%
        filter(
          ground_value == ground_key,
          !is.na(treatment_value),
          stringr::str_trim(treatment_value) != ""
        ) %>%
        count(treatment_value, name = "n_treatment") %>%
        arrange(desc(n_treatment), treatment_value)

      if (nrow(treatment_counts) == 0) {
        return(ground_row)
      }

      treatment_counts <- treatment_counts %>%
        filter(
          stringr::str_to_lower(stringr::str_trim(treatment_value)) !=
            stringr::str_to_lower(stringr::str_trim(ground_key))
        )

      if (nrow(treatment_counts) == 0) {
        return(ground_row)
      }

      treatment_rows <- treatment_counts %>%
        mutate(
          variable = paste0("Treatment group: ", treatment_value),
          n = format_n(n_treatment),
          mean = "",
          sd = "",
          min = "",
          median = "",
          max = "",
          percent = format_num2(100 * n_treatment / sum(n_treatment))
        ) %>%
        select(variable, n, mean, sd, min, median, max, percent)

      bind_rows(ground_row, treatment_rows)
    }
  )

  bind_rows(header_row, ground_rows)
}

categorical.map <- tibble::tribble(
  ~var_name, ~var_label,
  "region_agg", "Region (aggregated)",
  "info_binary", "Application Detail Level",
  "agent_type", "Type of Housing Provider",
  "matched_design", "Matched Design",
  "audit_type", "Mode of Application",
  "callback_type", "Callback Type",
  "airbnb", "Airbnb Indicator",
  "can_gender", "Candidate Gender",
  "can_education", "Candidate Education",
  "can_employment", "Candidate Employment",
  "can_migrant_generation", "Candidate Migrant Generation",
  "peer_reviewed", "Peer-reviewed Study",
  "language", "Study Language"
)

housing.meta.panel.b <- bind_rows(
  summarise_ground_treatment(housing.meta),
  purrr::pmap_dfr(
    categorical.map,
    function(var_name, var_label) summarise_categorical(housing.meta, var_name, var_label)
  )
)

panel.a.header <- tibble(
  variable = "Panel A. Continuous Variables",
  n = "",
  mean = "",
  sd = "",
  min = "",
  median = "",
  max = "",
  percent = ""
)

panel.b.header <- tibble(
  variable = "Panel B. Categorical Variables",
  n = "",
  mean = "",
  sd = "",
  min = "",
  median = "",
  max = "",
  percent = ""
)

housing.meta.table1 <- bind_rows(
  panel.a.header,
  housing.meta.panel.a,
  blank_row(),
  panel.b.header,
  housing.meta.panel.b
) %>%
  rename(
    Variable = variable,
    N = n,
    Mean = mean,
    SD = sd,
    Min = min,
    Median = median,
    Max = max,
    Percent = percent
  )

modelsummary::datasummary_df(
  housing.meta.table1,
  output = file.path(dir.descriptives, "meta_housing_table1_descriptives.html"),
  title = "Table 1. Descriptive statistics of focal variables"
)

# Appendix: Countries by region ####
housing.meta.ctr.map <- housing.meta %>%
  mutate(
    country_key = dplyr::coalesce(as.character(country_iso3), as.character(country)),
    country = as.character(country),
    region = as.character(region),
    region_agg = as.character(region_agg)
  ) %>%
  filter(
    !is.na(country_key), country_key != "",
    !is.na(region), region != "",
    !is.na(region_agg), region_agg != ""
  ) %>%
  distinct(country_key, country, region, region_agg)

housing.meta.ctr.map.check <- housing.meta.ctr.map %>%
  group_by(country_key) %>%
  summarise(
    n_region = n_distinct(region),
    n_region_agg = n_distinct(region_agg),
    .groups = "drop"
  ) %>%
  filter(n_region > 1 | n_region_agg > 1)

if (nrow(housing.meta.ctr.map.check) > 0) {
  stop(
    "Country-to-region mapping is not one-to-one for some countries. ",
    "Check country/region/region_agg coding in housing.meta before producing appendix table.",
    call. = FALSE
  )
}

housing.meta.ctr.stats <- housing.meta %>%
  mutate(
    country_key = dplyr::coalesce(as.character(country_iso3), as.character(country)),
    country = as.character(country),
    region = as.character(region),
    region_agg = as.character(region_agg)
  ) %>%
  filter(
    !is.na(country_key), country_key != "",
    !is.na(country), country != "",
    !is.na(region), region != "",
    !is.na(region_agg), region_agg != ""
  ) %>%
  group_by(country_key, country, region_agg, region) %>%
  summarise(k_effects = n(), .groups = "drop") %>%
  mutate(percent_effects = 100 * k_effects / sum(k_effects)) %>%
  arrange(region_agg, region, country)

housing.meta.appendix.countries <- purrr::map_dfr(
  unique(housing.meta.ctr.stats$region_agg),
  function(region_agg_value) {
    df_agg <- housing.meta.ctr.stats %>%
      filter(region_agg == region_agg_value)

    agg_header <- tibble(
      country = paste0("Region: ", region_agg_value),
      k = "",
      percent = ""
    )

    region_blocks <- purrr::map_dfr(
      unique(df_agg$region),
      function(region_value) {
        df_region <- df_agg %>%
          filter(region == region_value)

        country_rows <- df_region %>%
          transmute(
            country = country,
            k = format_n(k_effects),
            percent = format_num2(percent_effects)
          )

        country_rows
      }
    )

    bind_rows(
      agg_header,
      region_blocks,
      tibble(country = "", k = "", percent = "")
    )
  }
) %>%
  rename(
    Country = country,
    k = k,
    Percent = percent
  )

if (nrow(housing.meta.appendix.countries) > 0 && all(housing.meta.appendix.countries[nrow(housing.meta.appendix.countries), ] == "")) {
  housing.meta.appendix.countries <- housing.meta.appendix.countries[-nrow(housing.meta.appendix.countries), ]
}

modelsummary::datasummary_df(
  housing.meta.appendix.countries,
  output = file.path(dir.descriptives.appendix, "meta_housing_appendix_countries_by_region.html"),
  title = "Table A1. Number of effects by country"
)

# Appendix: Effects by period ####
housing.meta.appendix.period <- housing.meta %>%
  filter(!is.na(period_narrow)) %>%
  mutate(period_narrow = as.character(period_narrow)) %>%
  group_by(period_narrow) %>%
  summarise(
    k_effects = n(),
    k_studies = n_distinct(study),
    .groups = "drop"
  ) %>%
  arrange(period_narrow) %>%
  mutate(percent_effects = format_num2(100 * k_effects / sum(k_effects))) %>%
  transmute(
    Period  = period_narrow,
    k       = format_n(k_effects),
    Studies = format_n(k_studies),
    Percent = percent_effects
  )

modelsummary::datasummary_df(
  housing.meta.appendix.period,
  output = file.path(dir.descriptives.appendix, "meta_housing_appendix_period_k_effects.html"),
  title = "Table A. Number of effects and studies by 5-year period"
)

housing.meta.descriptives <- list(
  table1              = housing.meta.table1,
  countries_by_region = housing.meta.appendix.countries,
  period_k_effects    = housing.meta.appendix.period
)
saveRDS(housing.meta.descriptives, file = file.path(dir.descriptives, "housing_meta_descriptives.rds"))
writexl::write_xlsx(
  housing.meta.descriptives,
  path = file.path(dir.descriptives, "meta_housing_descriptives.xlsx")
)

message("Descriptive tables complete.")
