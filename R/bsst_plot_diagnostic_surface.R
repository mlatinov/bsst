## ============================================================
## bsst_plot_diagnostic_surface()
## Same 1D/2D layout, colored by sampler diagnostics instead of error
## ============================================================
#' Plot Sampler Diagnostics Across the Design Space
#'
#' Same 1D/2D layout as \code{\link{bsst_plot_stress_1d}} /
#' \code{\link{bsst_plot_stress_2d}}, but colored by convergence diagnostics
#' (percent divergent fits, or min ESS) instead of recovery error. Answers:
#' "is the region where recovery looks bad actually a region where the
#' sampler is failing, rather than a genuine bias/identifiability problem?"
#' — a distinct failure mode requiring a different fix (reparameterize,
#' longer warmup) than genuine parameter non-identifiability.
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var Character; design variable for the x-axis.
#' @param y_var Character or \code{NULL}; if supplied, plots a 2D heatmap;
#'   if \code{NULL}, plots a 1D line.
#' @param parameter Character; single tracked parameter name. Required
#'   (diagnostics are stored per-fit, but filtering to one parameter keeps
#'   the row structure consistent with the other plotting functions).
#' @param diagnostic_metric One of \code{"pct_divergent"}, \code{"rhat"},
#'   \code{"ess_bulk"}, \code{"ess_tail"}. \code{"pct_divergent"} is only
#'   meaningful when \code{x} has multiple rows per design cell (rare for
#'   \code{bsst_stress()}, which fits one replicate per point) — for a
#'   single-replicate design it collapses to a 0/100 indicator of whether
#'   that one fit had any divergences. Default \code{"rhat"}.
#' @param fix_others Named list of fixed values for other varying design
#'   variables (required if more than \code{x_var}/\code{y_var} vary).
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_diagnostic_surface <- function(x, x_var, y_var = NULL, parameter,
                                          diagnostic_metric = c("rhat", "ess_bulk", "ess_tail", "pct_divergent"),
                                          fix_others = NULL,
                                          palette = "dark_research",
                                          title = NULL) {

  diagnostic_metric <- match.arg(diagnostic_metric)
  is_2d <- !is.null(y_var)

  if (missing(parameter) || is.null(parameter) || length(parameter) != 1) {
    stop("`parameter` must be supplied explicitly as a single parameter name.")
  }
  if (!"raw_table" %in% names(x)) stop("x must be the output of bsst_stress() or bsst_bo().")

  df <- x$raw_table
  used_vars <- if (is_2d) c(x_var, y_var) else x_var
  dvars <- .bsst_get_design_vars(x)
  for (v in used_vars) if (!v %in% dvars$varying) stop("'", v, "' is not a varying design variable.")

  other_varying <- setdiff(dvars$varying, used_vars)
  if (length(other_varying) > 0) {
    missing_fix <- setdiff(other_varying, names(fix_others))
    if (length(missing_fix) > 0) {
      stop("The design has other varying variables besides ", paste(used_vars, collapse = ", "), ": ",
           paste(other_varying, collapse = ", "),
           ". You must supply `fix_others` covering: ", paste(missing_fix, collapse = ", "))
    }
  }

  df <- df[df$parameter == parameter, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows found for parameter '", parameter, "'.")

  if (length(other_varying) > 0) {
    for (nm in other_varying) {
      target <- fix_others[[nm]]
      if (is.numeric(df[[nm]])) {
        closest_val <- df[[nm]][which.min(abs(df[[nm]] - target))]
        if (closest_val != target) {
          message("fix_others: no exact match for '", nm, "' = ", target,
                  ". Using nearest available value: ", closest_val, ".")
        }
        df <- df[df[[nm]] == closest_val, , drop = FALSE]
      } else {
        matches <- df[[nm]] == target
        if (!any(matches)) stop("fix_others value '", target, "' for '", nm, "' has no matching rows.")
        df <- df[matches, , drop = FALSE]
      }
    }
  }
  if (nrow(df) == 0) stop("No rows remain after applying fix_others filters.")

  df$metric_val <- switch(diagnostic_metric,
    rhat = df$rhat,
    ess_bulk = df$ess_bulk,
    ess_tail = df$ess_tail,
    pct_divergent = ifelse(df$divergences > 0, 100, 0)
  )

  metric_label <- switch(diagnostic_metric,
    rhat = "Max R-hat", ess_bulk = "Min bulk ESS",
    ess_tail = "Min tail ESS", pct_divergent = "Fit had divergences (%)")

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette

  if (!is_2d) {
    df <- df[order(df[[x_var]]), ]
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = metric_val)) +
      ggplot2::geom_line(color = pal$text_secondary, linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(color = converged), size = 2.6) +
      ggplot2::scale_color_manual(values = c(`TRUE` = pal$ordinal[1], `FALSE` = pal$accent),
                                   name = "Converged") +
      ggplot2::labs(x = x_var, y = metric_label,
                    title = title %||% paste0(metric_label, " vs. ", x_var, " \u2014 ", parameter)) +
      ggtheme$theme
    return(p)
  }

  if (any(duplicated(df[, c(x_var, y_var)]))) {
    warning("Multiple rows share the same (", x_var, ", ", y_var, "); averaging ", diagnostic_metric, ".")
    df <- stats::aggregate(df$metric_val, by = list(xx = df[[x_var]], yy = df[[y_var]]), FUN = mean, na.rm = TRUE)
    names(df) <- c(x_var, y_var, "metric_val")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = metric_val)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = pal$sequential, name = metric_label) +
    ggplot2::labs(x = x_var, y = y_var,
                  title = title %||% paste0(metric_label, " \u2014 ", parameter)) +
    ggtheme$theme

  p
}