is_informative_fixed <- function(x) {
  if (is.numeric(x)) {
    return(sum(is.finite(x)) > 1)
  }
  dplyr::n_distinct(x[!is.na(x)]) > 1
}

term_df <- function(x) {
  if (is.numeric(x)) {
    return(1L)
  }
  as.integer(max(dplyr::n_distinct(x[!is.na(x)]) - 1L, 0L))
}

interaction_vars <- function(term) {
  term_clean <- gsub("\\s+", "", as.character(term))
  vars <- unlist(strsplit(term_clean, ":", fixed = TRUE), use.names = FALSE)
  unique(vars[vars != ""])
}

rhs_rank <- function(data, rhs_terms) {
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & rhs_terms != ""]
  rhs_text <- if (length(rhs_terms) == 0) {
    "1"
  } else {
    paste(c("1", rhs_terms), collapse = " + ")
  }

  mm <- tryCatch(
    stats::model.matrix(stats::as.formula(paste("~", rhs_text)), data = data),
    error = function(e) NULL
  )
  if (is.null(mm)) {
    return(NA_integer_)
  }
  as.integer(qr(mm)$rank)
}

select_interaction_terms <- function(data, interaction_terms, base_terms) {
  interaction_candidates <- unique(interaction_terms[!is.na(interaction_terms) & interaction_terms != ""])
  if (length(interaction_candidates) == 0) {
    return(character(0))
  }

  selected <- character(0)
  current_terms <- base_terms
  current_rank <- rhs_rank(data, current_terms)
  if (!is.finite(current_rank)) {
    return(character(0))
  }

  for (term in interaction_candidates) {
    vars <- interaction_vars(term)
    if (length(vars) < 2) {
      next
    }
    if (!all(vars %in% names(data))) {
      next
    }
    if (!all(vars %in% base_terms)) {
      next
    }
    if (!all(vapply(vars, function(v) is_informative_fixed(data[[v]]), logical(1)))) {
      next
    }

    new_rank <- rhs_rank(data, c(current_terms, term))
    if (is.finite(new_rank) && new_rank > current_rank) {
      selected <- c(selected, term)
      current_terms <- c(current_terms, term)
      current_rank <- new_rank
    }
  }

  selected
}

calc_design_df <- function(data, fixed_terms, random_terms,
                           interaction_terms = character(0),
                           random_df_method = c("count_all", "study_only")) {
  random_df_method <- match.arg(random_df_method)
  fixed_df <- sum(vapply(fixed_terms, function(v) term_df(data[[v]]), integer(1)))
  interaction_df <- 0L
  if (length(interaction_terms) > 0) {
    current_terms <- fixed_terms
    current_rank <- rhs_rank(data, current_terms)
    for (term in interaction_terms) {
      new_rank <- rhs_rank(data, c(current_terms, term))
      if (is.finite(current_rank) && is.finite(new_rank) && new_rank > current_rank) {
        interaction_df <- interaction_df + as.integer(new_rank - current_rank)
        current_terms <- c(current_terms, term)
        current_rank <- new_rank
      }
    }
  }
  random_df <- switch(random_df_method,
    count_all = length(random_terms),
    study_only = {
      study_levels <- dplyr::n_distinct(data$study[!is.na(data$study)])
      if ("study" %in% random_terms && study_levels > 1) 1L else 0L
    }
  )
  total_df <- 1L + fixed_df + interaction_df + random_df
  list(fixed_df = fixed_df, interaction_df = interaction_df, random_df = random_df, total_df = total_df)
}
