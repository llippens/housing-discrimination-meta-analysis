# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c(
  "here", "readxl", "readr", "janitor",
  "tidyr", "forcats", "stringr", "dplyr", "lubridate", "tibble"
)

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
dir.data <- file.path(dir.root, "1_data")
dir.public <- file.path(dir.root, "1_data", "public")
dir.tables <- file.path(dir.root, "4_tables")
dir.wrangling <- file.path(dir.tables, "1_wrangling")

ensure_dir(dir.public)
ensure_dir(dir.tables)
ensure_dir(dir.wrangling)

# Inputs ####
data.path <- file.path(dir.data, "260206_housing.xlsx")
housing.raw <- readxl::read_xlsx(data.path, sheet = "data") %>%
  janitor::clean_names()

required.cols <- c(
  "study", "cluster", "ground", "ground_abbr", "peer_reviewed",
  "treatment_group", "control_group", "can_gender", "can_education", "can_employment",
  "can_migrant_generation", "year_research_start", "quarter_research_start",
  "year_research_end", "quarter_research_end", "country", "country_iso3", "region",
  "subregion", "subregion2", "callback_type", "agent_type", "audit_type", "airbnb",
  "information", "matched_design", "language", "application_count_maj", "application_count_min",
  "callback_count_maj", "callback_count_min"
)
assert_required_cols(housing.raw, required.cols, "raw housing data")

ground.scope <- c("ero", "gen", "hed", "seo", "soc")
peer.check.cols <- names(housing.raw)[stringr::str_detect(names(housing.raw), "^peer_check")]
precomputed.effect.cols <- c(
  "response_rate",
  "positive_response_rate_maj",
  "positive_response_rate_min",
  "positive_response_ratio",
  "positive_response_ratio_ln",
  "positive_response_ratio_se_delta",
  "discrimination_ratio"
)

# Data prep ####
housing.data <- 
  housing.raw %>%
  mutate(
    cluster = as.character(cluster),
    study = as.character(study),
    study = case_when(
      !is.na(study) & study != "" ~ study,
      TRUE ~ paste0("row_", row_number())
    ),
    cluster = case_when(
      !is.na(cluster) & cluster != "" ~ cluster,
      TRUE ~ study
    ),
    ground = str_squish(as.character(ground)),
    ground = case_when(
      ground == "Race, ethnic identity, and national origin" ~ "Ethno-racial origin",
      ground == "Sex and gender" ~ "Gender",
      ground == "Sexual orientation" ~ "Sexual orientation",
      ground == "Wealth and class" ~ "Social origin",
      ground == "Health and disability" ~ "Health and disability",
      TRUE ~ ground
    ),
    ground_abbr = str_to_lower(str_squish(as.character(ground_abbr))),
    ground_abbr = case_when(
      ground_abbr == "ren" ~ "ero",
      ground_abbr == "seg" ~ "gen",
      ground_abbr == "hed" ~ "hed",
      ground_abbr == "seo" ~ "seo",
      ground_abbr == "wea" ~ "soc",
      TRUE ~ ground_abbr
    ),
    language = str_squish(as.character(language)),
    language = case_when(
      language %in% c("French", "German") ~ "Other",
      TRUE ~ language
    ),
    treatment_control = str_c(treatment_group, " x ", control_group),
    peer_reviewer_count = if (length(peer.check.cols) == 0) {
      0L
    } else {
      rowSums(
        across(
          all_of(peer.check.cols),
          ~ !is.na(.) & stringr::str_trim(as.character(.)) != ""
        )
      )
    }
  ) %>%
  relocate(treatment_control, .after = control_group) %>%
  mutate(
    peer_reviewed = factor(peer_reviewed, levels = c("No", "Yes")),
    matched_design = factor(matched_design, levels = c("No", "Yes")),
    callback_type = factor(callback_type, levels = c("Positive reaction", "Viewing invitation")),
    ground = factor(ground, levels = c("Ethno-racial origin", "Gender", "Health and disability", "Sexual orientation", "Social origin")),
    ground_abbr = factor(ground_abbr, levels = c("ero", "gen", "hed", "seo", "soc")),
    agent_type = case_when(
      agent_type %in% c("Mixed", "Unknown", "", NA_character_) ~ "Mixed or Unknown",
      TRUE ~ as.character(agent_type)
    ) %>% factor(levels = c("Private owner", "Real estate agent", "Mixed or Unknown")),
    audit_type = case_when(
      audit_type %in% c("Telephone", "Mixed") ~ "Telephone or Mixed",
      TRUE ~ as.character(audit_type)
    ) %>% factor(levels = c("Email", "Telephone or Mixed", "In person")),
    airbnb = factor(airbnb, levels = c("No", "Yes")),
    info_binary = case_when(
      information %in% c("No information", "", NA_character_) ~ "Standard",
      TRUE ~ "Extensive"
    ) %>% factor(levels = c("Standard", "Extensive")),
    information = as.factor(information),
    region = as.factor(region),
    region_agg = case_when(
      region %in% c("Europe", "North America") ~ as.character(region),
      TRUE ~ "Other"
    ) %>% factor(levels = c("Europe", "North America", "Other")),
    subregion = as.factor(subregion),
    subregion2 = as.factor(subregion2),
    country = {
      country_fac <- as.factor(country)
      if ("United States" %in% levels(country_fac)) {
        forcats::fct_relevel(country_fac, "United States")
      } else {
        country_fac
      }
    },
    country_iso3 = as.factor(country_iso3),
    language = {
      language_fac <- as.factor(language)
      if ("English" %in% levels(language_fac)) {
        forcats::fct_relevel(language_fac, "English")
      } else {
        language_fac
      }
    },
    can_gender = case_when(
      can_gender %in% c("Unknown", "", NA_character_) ~ "Unknown or NA",
      TRUE ~ can_gender
    ) %>% factor(levels = c("Male", "Female", "Both", "Unknown or NA")),
    can_education = case_when(
      can_education %in% c("Unknown", "", NA_character_, "No higher education", "Various") ~ "Various or Unknown",
      can_education %in% c("Higher education", "College", "Bachelor's", "University") ~ "Higher education",
      TRUE ~ can_education
    ) %>%
      factor() %>%
      forcats::fct_relevel("Higher education", "Various or Unknown"),
    can_employment = case_when(
      can_employment %in% c("", NA_character_, "Unknown", "Mixed") ~ "Mixed or Unknown",
      TRUE ~ can_employment
    ) %>%
      factor() %>%
      forcats::fct_relevel("Employed", "Unemployed", "Mixed or Unknown"),
    can_migrant_generation = case_when(
      ground_abbr != "ero" ~ "Unknown or NA",
      can_migrant_generation %in% c("", NA_character_, "Unknown") ~ "Unknown or NA",
      TRUE ~ can_migrant_generation
    ) %>%
      factor() %>%
      forcats::fct_relevel("First generation", "Second generation", "Unknown or NA"),
    treatment_group = if_else(
      ground_abbr == "ero",
      case_when(
        str_detect(str_to_lower(treatment_group), "asian") ~ "Asian",
        str_detect(str_to_lower(treatment_group), "europe|european|north\\s*american|northern\\s*american|white") ~ "Eastern European and other White",
        str_detect(str_to_lower(treatment_group), "arab|north\\s*african|northern\\s*african|turkish") ~ "Arab, North African, Turkish",
        str_detect(str_to_lower(treatment_group), "latin|hispanic") ~ "Hispanic, Southern American",
        str_detect(str_to_lower(treatment_group), "black") ~ "Black, Central African, African American",
        TRUE ~ "Other"
      ),
      treatment_group
    ) %>% as.factor()
  ) %>%
  mutate(
    quarter_research_start = suppressWarnings(as.numeric(quarter_research_start)),
    quarter_research_end = suppressWarnings(as.numeric(quarter_research_end)),
    year_research_start = suppressWarnings(as.numeric(year_research_start)),
    year_research_end = suppressWarnings(as.numeric(year_research_end)),
    time_research_start = if_else(
      is.na(quarter_research_start),
      year_research_start,
      year_research_start + (quarter_research_start - 1) / 4
    ),
    time_research_end = if_else(
      is.na(quarter_research_end),
      year_research_end,
      year_research_end + (quarter_research_end - 1) / 4
    ),
    year_mid = floor(time_research_start + (time_research_end - time_research_start) / 2),
    year_mid_z = as.numeric(scale(year_mid)),
    year_mid_fct = factor(year_mid),
    period_broad = case_when(
      year_mid >= 1999 & year_mid <= 2010 ~ "1999-2010",
      year_mid >= 2011 & year_mid <= 2020 ~ "2011-2020",
      year_mid >= 2021 & year_mid <= 2024 ~ "2021-2024",
      TRUE ~ NA_character_
    ) %>% factor(),
    period_narrow = case_when(
      year_mid >= 1999 & year_mid <= 2005 ~ "1999-2005",
      year_mid >= 2006 & year_mid <= 2010 ~ "2006-2010",
      year_mid >= 2011 & year_mid <= 2015 ~ "2011-2015",
      year_mid >= 2016 & year_mid <= 2020 ~ "2016-2020",
      year_mid >= 2021 & year_mid <= 2024 ~ "2021-2024",
      TRUE ~ NA_character_
    ) %>% factor()
  ) %>%
  select(-any_of(precomputed.effect.cols)) %>%
  mutate(
    noncallback_count_maj = application_count_maj - callback_count_maj,
    noncallback_count_min = application_count_min - callback_count_min,
    rrate = (callback_count_maj + callback_count_min) / (application_count_maj + application_count_min),
    rrate_ln = log(rrate),
    prrate_maj = callback_count_maj / application_count_maj,
    prrate_min = callback_count_min / application_count_min,
    prratio = prrate_min / prrate_maj,
    prratio_ln = log(prratio),
    prratio_ln_se = sqrt(
      (1 / callback_count_maj) +
        (1 / callback_count_min) -
        (1 / application_count_maj) -
        (1 / application_count_min)
    ),
    prratio_ln_v = prratio_ln_se^2,
    prratio_ln_p = 1 / prratio_ln_se,
    prratio_ln_v_inv = 1 / prratio_ln_v,
    prratio_std = prratio_ln / prratio_ln_se,
    dratio = prratio - 1,
    odds_maj = callback_count_maj / noncallback_count_maj,
    odds_min = callback_count_min / noncallback_count_min,
    oddsratio = odds_min / odds_maj,
    oddsratio_ln = log(oddsratio),
    oddsratio_ln_se = sqrt(
      (1 / callback_count_maj) +
        (1 / callback_count_min) +
        (1 / (application_count_maj - callback_count_maj)) +
        (1 / (application_count_min - callback_count_min))
    )
  ) %>%
  mutate(
    flag_in_scope_ground = ground_abbr %in% ground.scope,
    flag_complete_counts = !is.na(application_count_maj) &
      !is.na(application_count_min) &
      !is.na(callback_count_maj) &
      !is.na(callback_count_min),
    flag_valid_effect = is.finite(prratio_ln) & is.finite(prratio_ln_v) & prratio_ln_v > 0,
    flag_main_sample = flag_in_scope_ground &
      flag_complete_counts &
      flag_valid_effect &
      audit_type != "In person"
  )

# Guardrails: Ground taxonomy consistency ####
ground.map.expected <- tibble::tribble(
  ~ground, ~ground_abbr,
  "Ethno-racial origin", "ero",
  "Gender", "gen",
  "Health and disability", "hed",
  "Sexual orientation", "seo",
  "Social origin", "soc"
)

ground.map.observed <- housing.data %>%
  mutate(
    ground = as.character(ground),
    ground_abbr = as.character(ground_abbr)
  ) %>%
  filter(ground_abbr %in% ground.scope) %>%
  distinct(ground, ground_abbr)

ground.map.unexpected <- ground.map.observed %>%
  anti_join(ground.map.expected, by = c("ground", "ground_abbr"))

if (nrow(ground.map.unexpected) > 0) {
  stop(
    "Unexpected ground mapping found in in-scope data: ",
    paste0(
      ground.map.unexpected$ground, " -> ", ground.map.unexpected$ground_abbr,
      collapse = "; "
    ),
    call. = FALSE
  )
}

# Analysis samples ####
housing.meta.with.inperson <- housing.data %>%
  filter(flag_in_scope_ground, flag_complete_counts, flag_valid_effect) %>%
  droplevels()

housing.meta <- housing.meta.with.inperson %>%
  filter(audit_type != "In person") %>%
  select(-any_of(peer.check.cols), -any_of(c("issue_resolved", "can_age_mid"))) %>%
  relocate(study, cluster, .before = 1) %>%
  relocate(treatment_control, .after = control_group) %>%
  relocate(
    noncallback_count_maj, noncallback_count_min,
    rrate, rrate_ln,
    prrate_maj, prrate_min, prratio, prratio_ln, prratio_ln_se,
    prratio_ln_v, prratio_ln_p, prratio_ln_v_inv, prratio_std, dratio,
    odds_maj, odds_min, oddsratio, oddsratio_ln, oddsratio_ln_se,
    .after = callback_count_min
  ) %>%
  droplevels()

# Coverage summaries ####
housing.meta.inclusion.flow <-
  tibble::tibble(
    step = c(
      "Raw rows",
      "In manuscript ground scope",
      "Complete core counts",
      "Valid effect-size rows",
      "Main sample (exclude in-person)"
    ),
    rows = c(
      nrow(housing.data),
      sum(housing.data$flag_in_scope_ground, na.rm = TRUE),
      sum(housing.data$flag_in_scope_ground & housing.data$flag_complete_counts, na.rm = TRUE),
      sum(housing.data$flag_in_scope_ground & housing.data$flag_complete_counts & housing.data$flag_valid_effect, na.rm = TRUE),
      sum(housing.data$flag_main_sample, na.rm = TRUE)
    ),
    studies = c(
      n_distinct(housing.data$study),
      n_distinct(housing.data$study[housing.data$flag_in_scope_ground]),
      n_distinct(housing.data$study[housing.data$flag_in_scope_ground & housing.data$flag_complete_counts]),
      n_distinct(housing.data$study[housing.data$flag_in_scope_ground & housing.data$flag_complete_counts & housing.data$flag_valid_effect]),
      n_distinct(housing.data$study[housing.data$flag_main_sample])
    )
  )

housing.meta.ground.coverage <-
  housing.meta %>%
  group_by(ground, ground_abbr) %>%
  summarise(
    rows = n(),
    studies = n_distinct(study),
    min_year = min(year_mid, na.rm = TRUE),
    max_year = max(year_mid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(rows))

# Outputs ####
readr::write_csv(housing.meta.inclusion.flow, file.path(dir.wrangling, "meta_housing_inclusion_flow.csv"), na = "")
readr::write_csv(housing.meta.ground.coverage, file.path(dir.wrangling, "meta_housing_ground_coverage.csv"), na = "")
saveRDS(housing.meta, file.path(dir.public, "housing_meta.rds"))
readr::write_csv(housing.meta, file.path(dir.public, "housing_meta.csv"), na = "")

message("Wrangling complete: ", nrow(housing.meta), " analysable effect rows across ", n_distinct(housing.meta$study), " studies.")
