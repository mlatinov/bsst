## ============================================================
## bsst_plot_stress_1d()
## ============================================================
#' Plot Recovery/Stress Metric Against One Design Variable
#'
#' Visualizes how a recovery metric (error, standardized error, or stress
#' score) changes across one simulation design variable, for a single
#' tracked parameter. Answers: "as I vary this simulation setting, at what
#' point does recovery for this parameter stabilize or break down?"
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var Character; name of the design variable to place on the
#'   x-axis (must be a varying variable in \code{x}).
#' @param parameter Character; name of a single tracked parameter to plot.
#'   Required — must be supplied explicitly (no default / auto-facet).
#' @param y_var One of \code{"std_error"}, \code{"error"},
#'   \code{"stress_score"}. Default \code{"std_error"}.
#' @param fix_others Named list giving fixed values for any other varying
#'   design variables present in \code{x}, required whenever more than
#'   \code{x_var} varies in the underlying design. Rows not matching these
#'   fixed values exactly are excluded.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param show_nonconverged Logical; if \code{TRUE} (default), points that
#'   failed convergence diagnostics are still shown but marked with a
#'   distinct shape, rather than silently dropped — since a pattern of
#'   non-convergence across the design space is itself diagnostic
#'   information.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_stress_1d <- function(x, x_var, parameter,
                                 y_var = c("std_error", "error", "stress_score"),
                                 fix_others = NULL,
                                 palette = "dark_research",
                                 show_nonconverged = TRUE,
                                 title = NULL) {

  y_var <- match.arg(y_var)

  if (missing(parameter) || is.null(parameter)) {
    stop("`parameter` must be supplied explicitly (no default / auto-facet). ",
         "Pass a single tracked parameter name.")
  }
  if (length(parameter) != 1) {
    stop("`parameter` must be a single parameter name for bsst_plot_stress_1d(). ",
         "Call this function once per parameter if you need multiple.")
  }
  if (!"raw_table" %in% names(x)) {
    stop("x must be the output of bsst_stress() or bsst_bo().")
  }

  df <- x$raw_table
  if (!x_var %in% names(df)) {
    stop("x_var '", x_var, "' not found in x$raw_table.")
  }
  if (!y_var %in% names(df)) {
    stop("y_var '", y_var, "' not found in x$raw_table (was objective_fn NULL when this was generated?).")
  }

  dvars <- .bsst_get_design_vars(x)
  if (!x_var %in% dvars$varying) {
    stop("x_var '", x_var, "' is not a varying design variable in x.")
  }

  other_varying <- setdiff(dvars$varying, x_var)

  if (length(other_varying) > 0) {
    missing_fix <- setdiff(other_varying, names(fix_others))
    if (length(missing_fix) > 0) {
      stop("The design has other varying variables besides '", x_var, "': ",
           paste(other_varying, collapse = ", "),
           ". You must supply `fix_others` covering: ",
           paste(missing_fix, collapse = ", "),
           " (e.g. fix_others = list(", missing_fix[1], " = <value>)).")
    }
  }

  ## ---- filter to the requested parameter ----
  df <- df[df$parameter == parameter, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No rows found for parameter '", parameter, "' in x$raw_table.")
  }

  ## ---- filter to fixed values of other varying variables ----
  if (length(other_varying) > 0) {
    for (nm in other_varying) {
      target <- fix_others[[nm]]
      matches <- df[[nm]] == target
      if (!any(matches)) {
        available <- paste(sort(unique(df[[nm]])), collapse = ", ")
        stop("fix_others value ", target, " for '", nm, "' has no matching rows. ",
             "Available values: ", available)
      }
      df <- df[matches, , drop = FALSE]
    }
  }

  if (nrow(df) == 0) {
    stop("No rows remain after applying fix_others filters.")
  }

  th <- .bsst_get_palette(palette)
  ggtheme <- .bsst_theme_ggplot(palette)

  df <- df[order(df[[x_var]]), ]
  df$converged_label <- ifelse(df$converged, "Converged", "Diagnostics failed")

  y_label <- switch(y_var,
    std_error    = "Standardized error  (mean\u2212true) / posterior sd",
    error        = "Error  (mean\u2212true)",
    stress_score = "Stress score"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]))

  if (y_var != "stress_score") {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = th$grid, linewidth = 0.6, linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_line(color = th$accent, linewidth = 0.8, alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(shape = converged_label),
                         color = th$accent, size = 2.6) +
    ggplot2::scale_shape_manual(values = c("Converged" = 16, "Diagnostics failed" = 4),
                                 name = NULL, drop = show_nonconverged) +
    ggplot2::labs(
      x = x_var, y = y_label,
      title = title %||% paste0(y_label, " vs. ", x_var, " \u2014 ", parameter),
      subtitle = if (length(other_varying) > 0) {
        paste0("Fixed: ", paste(sprintf("%s = %s", names(fix_others), unlist(fix_others)), collapse = ", "))
      } else NULL
    ) +
    ggtheme$theme

  if (!show_nonconverged) {
    p$data <- p$data[p$data$converged, ]
  }

  p
}