## ============================================================
## bsst_plot_profile()
## Thin wrapper over bsst_plot_stress_1d() — "vary one, fix/average others"
## ============================================================
#' Plot a Recovery Profile Across One Design Variable
#'
#' A convenience wrapper around \code{\link{bsst_plot_stress_1d}} for the
#' common "profile" use case: vary one simulation design variable while
#' holding all others fixed at specified values (or, when requested,
#' averaging over them instead of requiring an exact fix). Answers the same
#' question as \code{bsst_plot_stress_1d()} — "as I vary this setting, where
#' does recovery break down?" — with an explicit averaging option for users
#' who don't want to commit to one fixed value of the other design
#' variables.
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var Character; the design variable to place on the x-axis.
#' @param parameter Character; single tracked parameter name. Required.
#' @param y_var One of \code{"std_error"}, \code{"error"}, \code{"stress_score"}.
#' @param fix_others Named list of fixed values for other varying design
#'   variables. Required unless \code{average_others = TRUE}.
#' @param average_others Logical. If \code{TRUE}, instead of requiring
#'   \code{fix_others}, the function averages \code{y_var} across all other
#'   varying design variables at each level of \code{x_var} (via
#'   \code{mean()}, ignoring \code{NA}). This trades precision for not
#'   requiring the user to commit to specific values of nuisance design
#'   variables — useful for a first look, but note it can mask
#'   interactions that \code{bsst_plot_stress_2d()} would reveal (e.g., if
#'   the model is fine at low measurement error and badly biased at high
#'   measurement error, averaging over that dimension will wash the effect
#'   out into a misleadingly moderate line). Default \code{FALSE}.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param show_nonconverged Logical; see \code{\link{bsst_plot_stress_1d}}.
#'   Ignored (points always shown) when \code{average_others = TRUE}, since
#'   convergence status can't be meaningfully averaged.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_profile <- function(x, x_var, parameter,
                               y_var = c("std_error", "error", "stress_score"),
                               fix_others = NULL,
                               average_others = FALSE,
                               palette = "dark_research",
                               show_nonconverged = TRUE,
                               title = NULL) {

  y_var <- match.arg(y_var)

  if (!average_others) {
    ## delegate entirely to bsst_plot_stress_1d() — same semantics, same errors
    return(bsst_plot_stress_1d(
      x = x, x_var = x_var, parameter = parameter, y_var = y_var,
      fix_others = fix_others, palette = palette,
      show_nonconverged = show_nonconverged, title = title
    ))
  }

  ## ---- average_others = TRUE path ----
  if (missing(parameter) || is.null(parameter) || length(parameter) != 1) {
    stop("`parameter` must be supplied explicitly as a single parameter name.")
  }
  if (!"raw_table" %in% names(x)) {
    stop("x must be the output of bsst_stress() or bsst_bo().")
  }

  df <- x$raw_table
  if (!x_var %in% names(df)) stop("x_var '", x_var, "' not found in x$raw_table.")
  if (!y_var %in% names(df)) stop("y_var '", y_var, "' not found in x$raw_table.")

  dvars <- .bsst_get_design_vars(x)
  if (!x_var %in% dvars$varying) stop("x_var '", x_var, "' is not a varying design variable.")

  df <- df[df$parameter == parameter, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows found for parameter '", parameter, "'.")

  agg <- stats::aggregate(
    df[[y_var]], by = list(x_val = df[[x_var]]), FUN = mean, na.rm = TRUE
  )
  names(agg) <- c(x_var, y_var)
  agg <- agg[order(agg[[x_var]]), ]

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette

  y_label <- switch(y_var,
    std_error = "Standardized error  (mean\u2212true) / posterior sd (averaged)",
    error = "Error  (mean\u2212true) (averaged)",
    stress_score = "Stress score (averaged)")

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]))

  if (y_var != "stress_score") {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = pal$grid, linewidth = 0.6, linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_line(color = pal$accent, linewidth = 0.8) +
    ggplot2::geom_point(color = pal$accent, size = 2.6) +
    ggplot2::labs(
      x = x_var, y = y_label,
      title = title %||% paste0(y_label, " vs. ", x_var, " \u2014 ", parameter),
      subtitle = "Averaged across all other varying design variables \u2014 may mask interactions; see bsst_plot_stress_2d() to check"
    ) +
    ggtheme$theme

  p
}
