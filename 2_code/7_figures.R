# Packages ####
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

pkgs <- c("here", "dplyr", "ggplot2", "tibble", "forcats", "colorspace", "patchwork",
          "scales")
missing.pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing.pkgs) > 0) {
  pak::pkg_install(missing.pkgs)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# Font loading ####
if (!requireNamespace("showtext", quietly = TRUE)) pak::pkg_install("showtext")
library(showtext)
font_add_google("Atkinson Hyperlegible", family = "Atkinson Hyperlegible Next")
showtext_auto()

# Helpers ####
source(file.path(here::here(), "2_code", "functions", "core", "f_check_cols.R"))

# Reproducibility ####
set.seed(8888)

# Aesthetics ####
ugblue      <- "#1E64C8"
ugblue.dark <- "#123E7A"
ugred       <- "#DC4E28"
family      <- "Atkinson Hyperlegible Next"

est.colours <- c("RE" = ugblue, "UWLS" = "grey60")
est.fills   <- c("RE" = ugblue, "UWLS" = "grey60")

common_theme <-
  theme_minimal() +
  theme(text = element_text(family = family, size = 14, colour = "black"),
        panel.grid.major = element_line(linewidth = .25, colour = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 5, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 5, unit = "pt"),
                                    angle = 0, hjust = 0),
        axis.ticks = element_line(linewidth = .5, colour = "black"),
        panel.spacing = unit(5, units = "pt"),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank())

# Directories ####
dir.root      <- here::here()
dir.tables    <- file.path(dir.root, "4_tables")
dir.figures   <- file.path(dir.root, "3_figures")
dir.meta.main <- file.path(dir.tables, "2_meta_analysis")

dir.figures.appendix <- file.path(dir.figures, "appendix")

ensure_dir(dir.tables)
ensure_dir(dir.figures)
ensure_dir(dir.figures.appendix)

# Inputs ####
main.bundle.path <- file.path(dir.meta.main, "housing_meta_main.rds")

housing.meta.main.bundle <- {
  if (!file.exists(main.bundle.path))
    stop("Missing meta-analysis bundle: run 2_code/2_meta_analysis.R first.", call. = FALSE)
  readRDS(main.bundle.path)
}

housing.meta.table.main <- {
  tbl <- housing.meta.main.bundle[["housing.meta.summary"]]
  if (is.null(tbl))
    stop("Missing RE ground summary in bundle: run 2_code/2_meta_analysis.R first.", call. = FALSE)
  tbl
}

housing.meta.table.tgroup <- {
  tbl <- housing.meta.main.bundle[["housing.meta.treatment.summary"]]
  if (is.null(tbl))
    stop("Missing RE treatment-group summary in bundle: run 2_code/2_meta_analysis.R first.", call. = FALSE)
  tbl
}

housing.meta.table.pet.peese <- {
  tbl <- housing.meta.main.bundle[["housing.meta.pet.peese.selected"]]
  if (is.null(tbl))
    stop("Missing PET/PEESE selected table in bundle: run 2_code/2_meta_analysis.R first.", call. = FALSE)
  tbl
}

# Figure 2 — Ground-level (RE + UWLS-MRA) ####
ground.uwls <- housing.meta.table.pet.peese %>%
  filter(cell_type == "ground", !is.na(estimate_rr)) %>%
  transmute(ground = as.character(ground), estimator = "UWLS",
            est = estimate_rr, lo = conf_low_rr, hi = conf_high_rr)

ground.re <- housing.meta.table.main %>%
  filter(!is.na(rr_mv_cr2)) %>%
  transmute(ground = as.character(ground), estimator = "RE",
            est = rr_mv_cr2, lo = rr_mv_cr2_low, hi = rr_mv_cr2_high)

ground.plot <- bind_rows(ground.uwls, ground.re) %>%
  mutate(
    estimator = factor(estimator, levels = c("RE", "UWLS")),
    ground = fct_reorder(ground, ifelse(estimator == "RE", est, NA_real_),
                         .fun = function(x) mean(x, na.rm = TRUE), .na_rm = TRUE)
  )

p.ground <- ggplot(ground.plot, aes(x = est, y = ground,
                                     colour = estimator, fill = estimator)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey75", linewidth = 0.25) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.15, linewidth = 0.6,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, shape = 23, colour = "white", stroke = 0.4,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = est.colours) +
  scale_fill_manual(values = est.fills) +
  scale_x_continuous("Discrimination ratio") +
  scale_y_discrete("") +
  common_theme

# Figure 3 — Treatment-group level (RE + UWLS-MRA), panelled by ground ####
tgroup.uwls <- housing.meta.table.pet.peese %>%
  filter(cell_type == "treatment_group", !is.na(estimate_rr)) %>%
  transmute(ground = as.character(ground), ground_abbr = as.character(ground_abbr),
            treatment_group = as.character(treatment_group),
            estimator = "UWLS",
            est = estimate_rr, lo = conf_low_rr, hi = conf_high_rr)

tgroup.re <- housing.meta.table.tgroup %>%
  filter(!is.na(rr_mv_cr2)) %>%
  transmute(ground = as.character(ground), ground_abbr = as.character(ground_abbr),
            treatment_group = as.character(treatment_group),
            estimator = "RE",
            est = rr_mv_cr2, lo = rr_mv_cr2_low, hi = rr_mv_cr2_high)

tgroup.re.order <- tgroup.re %>%
  arrange(ground, est) %>%
  pull(treatment_group)

tgroup.plot <- bind_rows(tgroup.uwls, tgroup.re) %>%
  mutate(
    estimator = factor(estimator, levels = c("RE", "UWLS")),
    treatment_group = factor(treatment_group, levels = unique(tgroup.re.order))
  )

# Shared x limits across both figures ####
shared_xlim <- range(
  c(ground.plot$lo, ground.plot$hi, tgroup.plot$lo, tgroup.plot$hi),
  na.rm = TRUE
)

# Apply shared x limits to p.ground and save standalone ####
p.ground <- p.ground +
  coord_cartesian(xlim = shared_xlim) +
  guides(colour = guide_legend(override.aes = list(shape = NA, linetype = "solid", linewidth = 0.8)),
         fill = "none")

ggplot2::ggsave(
  filename = "meta_housing_ground_forest.png",
  path = dir.figures,
  plot = p.ground,
  width = 16,
  height = 9,
  units = "cm",
  dpi = 500
)

p.tgroup <- ggplot(tgroup.plot, aes(x = est, y = treatment_group,
                                     colour = estimator, fill = estimator)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey75", linewidth = 0.25) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.1, linewidth = 0.6,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, shape = 23, colour = "white", stroke = 0.4,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = est.colours) +
  scale_fill_manual(values = est.fills) +
  scale_x_continuous("Discrimination ratio") +
  scale_y_discrete("") +
  facet_grid(ground ~ ., scales = "free_y", space = "free_y") +
  common_theme +
  theme(strip.text = element_blank(),
        strip.background = element_blank()) +
  coord_cartesian(xlim = shared_xlim) +
  guides(colour = guide_legend(override.aes = list(shape = NA, linetype = "solid", linewidth = 0.8)),
         fill = "none")

ggplot2::ggsave(
  filename = "meta_housing_treatment_forest.png",
  path = dir.figures,
  plot = p.tgroup,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 500
)

# Figure — Combined panels (A: ground, B: treatment group) ####
n_ground_levels  <- nlevels(droplevels(ground.plot$ground))
n_tgroup_levels  <- nlevels(droplevels(tgroup.plot$treatment_group))

p.panels <- (p.ground / p.tgroup) +
  plot_layout(heights = c(n_ground_levels, n_tgroup_levels)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = family, size = 14, face = "bold", colour = "black"),
        legend.position = "top")

ggplot2::ggsave(
  filename = "meta_housing_main_panels.png",
  path = dir.figures,
  plot = p.panels,
  width = 18,
  height = 20,
  units = "cm",
  dpi = 500
)

# Funnel plots ####
z90 <- qnorm(.95)
z95 <- qnorm(.975)
z99 <- qnorm(.995)

filter_valid_effects <- function(data) {
  data %>%
    dplyr::filter(
      !is.na(prratio_ln_v),
      is.finite(prratio_ln_v),
      prratio_ln_v > 0,
      !is.na(prratio_ln),
      is.finite(prratio_ln)
    ) %>%
    dplyr::mutate(prratio_ln_se = sqrt(prratio_ln_v))
}

make_funnel_plot <- function(data_valid, pooled_ln, plot_label) {
  sei    <- data_valid$prratio_ln_se
  se_max <- ceiling(max(sei) * 10) / 10

  triangle90p <- tibble::tibble(
    group     = c(1, 1, 1),
    polygon.x = c(exp(0 + z90 * se_max), exp(0 - z90 * se_max), exp(0)),
    polygon.y = c(se_max, se_max, 0)
  )
  triangle95p <- tibble::tibble(
    group     = c(1, 1, 1),
    polygon.x = c(exp(0 + z95 * se_max), exp(0 - z95 * se_max), exp(0)),
    polygon.y = c(se_max, se_max, 0)
  )
  triangle99p <- tibble::tibble(
    group     = c(1, 1, 1),
    polygon.x = c(exp(0 + z99 * se_max), exp(0 - z99 * se_max), exp(0)),
    polygon.y = c(se_max, se_max, 0)
  )

  line95l <- tibble::tibble(
    x = c(exp(pooled_ln), exp(pooled_ln - z95 * se_max)), y = c(0, se_max)
  )
  line95r <- tibble::tibble(
    x = c(exp(pooled_ln), exp(pooled_ln + z95 * se_max)), y = c(0, se_max)
  )
  line99l <- tibble::tibble(
    x = c(exp(pooled_ln), exp(pooled_ln - z99 * se_max)), y = c(0, se_max)
  )
  line99r <- tibble::tibble(
    x = c(exp(pooled_ln), exp(pooled_ln + z99 * se_max)), y = c(0, se_max)
  )
  linese <- tibble::tibble(
    x = c(exp(pooled_ln), exp(pooled_ln)), y = c(0, se_max)
  )

  pts <- data_valid %>%
    dplyr::mutate(
      te   = exp(prratio_ln),
      sete = prratio_ln_se,
      w    = 1 / prratio_ln_v
    )

  log_lo   <- min(log(pts$te), pooled_ln) - 0.1
  log_hi   <- max(log(pts$te), pooled_ln) + 0.1
  x_breaks <- exp(pretty(c(log_lo, log_hi), n = 5))

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = triangle99p,
      mapping = ggplot2::aes(x = polygon.x, y = polygon.y, group = group),
      fill = "gray70"
    ) +
    ggplot2::geom_polygon(
      data = triangle95p,
      mapping = ggplot2::aes(x = polygon.x, y = polygon.y, group = group),
      fill = "gray85"
    ) +
    ggplot2::geom_polygon(
      data = triangle90p,
      mapping = ggplot2::aes(x = polygon.x, y = polygon.y, group = group),
      fill = "white"
    ) +
    ggplot2::geom_point(
      data = pts,
      mapping = ggplot2::aes(x = te, y = sete, size = w),
      alpha = 0.5, colour = "gray20", fill = "gray20"
    ) +
    ggplot2::scale_size_continuous(range = c(0.75, 3)) +
    ggplot2::geom_line(
      data = line95l,
      mapping = ggplot2::aes(x = x, y = y),
      colour = "gray20", linetype = "dotted"
    ) +
    ggplot2::geom_line(
      data = line95r,
      mapping = ggplot2::aes(x = x, y = y),
      colour = "gray20", linetype = "dotted"
    ) +
    ggplot2::geom_line(
      data = linese,
      mapping = ggplot2::aes(x = x, y = y),
      colour = "gray20", linetype = "dotted"
    ) +
    ggplot2::geom_line(
      data = line99l,
      mapping = ggplot2::aes(x = x, y = y),
      colour = "gray20", linetype = "dashed"
    ) +
    ggplot2::geom_line(
      data = line99r,
      mapping = ggplot2::aes(x = x, y = y),
      colour = "gray20", linetype = "dashed"
    ) +
    ggplot2::scale_x_log10(
      breaks = x_breaks,
      labels = scales::number_format(accuracy = 0.01)
    ) +
    ggplot2::scale_y_reverse(
      breaks = pretty(c(0, se_max)),
      labels = scales::comma_format(accuracy = 0.01)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(exp(log_lo), exp(log_hi)),
      ylim = c(se_max, 0)
    ) +
    ggplot2::labs(subtitle = plot_label, x = "Positive response ratio", y = "Standard error") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(family = family, size = 14, colour = "black"),
      panel.background = ggplot2::element_rect(fill = "gray95", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 8, r = 0, b = 0, l = 0), colour = "gray20"
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(t = 0, r = 8, b = 0, l = 0),
        angle = 90, hjust = 0.5, colour = "gray20"
      ),
      plot.subtitle = ggplot2::element_text(size = 11, colour = "gray20", hjust = 0.5),
      legend.position = "none"
    )
}

# Build funnel subsets ####
ground.pooled.ln <- setNames(
  log(housing.meta.table.main$rr_mv_cr2),
  housing.meta.table.main$ground_abbr
)

funnel.subs <- purrr::imap(
  housing.meta.main.bundle[["housing.meta.models"]],
  function(obj, g) {
    gvalid <- filter_valid_effects(obj$data)
    list(
      data_valid = gvalid,
      pooled_ln  = if (g %in% names(ground.pooled.ln) && is.finite(ground.pooled.ln[[g]]))
                     ground.pooled.ln[[g]] else
                     weighted.mean(gvalid$prratio_ln, 1 / gvalid$prratio_ln_v),
      label      = as.character(obj$data$ground[[1]])
    )
  }
)

funnel.tgroup.subs <- purrr::imap(
  housing.meta.main.bundle[["housing.meta.models"]],
  function(obj, g) {
    tgs <- unique(as.character(obj$data$treatment_group))
    tgs <- tgs[!is.na(tgs)]
    purrr::map(stats::setNames(tgs, tgs), function(tg) {
      tg_data  <- dplyr::filter(obj$data, as.character(treatment_group) == tg)
      tg_valid <- filter_valid_effects(tg_data)
      tg_row   <- dplyr::filter(
        housing.meta.table.tgroup,
        as.character(ground_abbr) == g,
        as.character(treatment_group) == tg
      )
      pooled_ln <- if (nrow(tg_row) > 0 && !is.na(tg_row$rr_mv_cr2[[1]]) &&
                       tg_row$rr_mv_cr2[[1]] > 0) {
        log(tg_row$rr_mv_cr2[[1]])
      } else if (nrow(tg_valid) > 0) {
        weighted.mean(tg_valid$prratio_ln, 1 / tg_valid$prratio_ln_v)
      } else {
        NA_real_
      }
      list(
        data_valid  = tg_valid,
        pooled_ln   = pooled_ln,
        label       = tg,
        ground_abbr = g,
        tg_safe     = tolower(gsub("[^a-zA-Z0-9]+", "_", tg))
      )
    })
  }
)
# Save ground funnel plots ####
for (sub_name in names(funnel.subs)) {
  sub <- funnel.subs[[sub_name]]
  if (nrow(sub$data_valid) < 3) {
    message("Skipping ground funnel: ", sub_name, " (< 3 valid effects)"); next
  }
  plt <- make_funnel_plot(sub$data_valid, sub$pooled_ln, sub$label)
  ggplot2::ggsave(
    filename = paste0("meta_housing_funnel_", sub_name, ".png"),
    path = dir.figures.appendix, plot = plt, device = "png",
    width = 15, height = 12, units = "cm", dpi = 1000, bg = "white"
  )
}

# Save treatment-group funnel plots (panelled by ground) ####
for (g in names(funnel.tgroup.subs)) {
  cells       <- funnel.tgroup.subs[[g]]
  valid_cells <- purrr::keep(cells, function(cell) {
    nrow(cell$data_valid) >= 3 && !is.null(cell$pooled_ln) && is.finite(cell$pooled_ln)
  })
  if (length(valid_cells) <= 1L) {
    message("Skipping tgroup panel: ", g, " (fewer than 2 qualifying cells, ground funnel sufficient)"); next
  }
  plts     <- purrr::map(valid_cells, function(cell) {
    make_funnel_plot(cell$data_valid, cell$pooled_ln, cell$label)
  })
  n_row    <- ceiling(length(plts) / 2L)
  combined <- patchwork::wrap_plots(plts, ncol = 2)
  ggplot2::ggsave(
    filename = paste0("meta_housing_funnel_tgroup_", g, ".png"),
    path = dir.figures.appendix, plot = combined, device = "png",
    width = 30, height = n_row * 12, units = "cm", dpi = 1000, bg = "white"
  )
}

message("Figures exported.")
