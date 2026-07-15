## ============================================================
## bsst_plot_recovery()
## Forest-plot: nested CI bars per parameter + true value marker
## ============================================================
#' Plot Parameter Recovery Intervals
#'
#' Visualizes the output of \code{\link{bsst_recover}} as a forest plot:
#' one row per tracked parameter, showing the posterior credible interval
#' (at the classification-relevant width) with the true value marked,
#' colored by recovery status.
#'
#' @param x Output list from \code{\link{bsst_recover}}.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param show_all_ci_levels Logical; if \code{TRUE}, draws all nested
#'   interval widths (thin-to-thick bars) instead of just the interval that
#'   determined classification. Default \code{FALSE}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_recovery <- function(x, palette = "dark_research",
                                show_all_ci_levels = FALSE, title = NULL) {

  if (!"summary_table" %in% names(x)) {
    stop("x must be the output of bsst_recover().")
  }
  df <- x$summary_table
  th <- .bsst_theme_ggplot(palette)

  df$recovery_status <- factor(df$recovery_status, levels = .bsst_status_levels)
  status_colors <- setNames(th$palette$ordinal, .bsst_status_levels)

  # order rows by severity for readability (worst at top)
  df$parameter <- factor(df$parameter, levels = rev(df$parameter[order(abs(df$std_error))]))

  p <- ggplot2::ggplot(df, ggplot2::aes(y = parameter))

  if (show_all_ci_levels) {
    # requires the raw fit's draws to recompute all levels; if only the
    # narrowest is stored, fall back to it with a warning
    warning("show_all_ci_levels = TRUE requires re-deriving nested intervals; ",
            "only the classification-relevant interval is stored in summary_table. ",
            "Falling back to single-interval display.")
  }

  p <- p +
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lower, xend = ci_upper, yend = parameter,
                   color = recovery_status),
      linewidth = 3, alpha = 0.85
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = estimate_mean), color = th$palette$text_main,
      size = 2.2, shape = 18
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = true_value), color = th$palette$accent,
      size = 3, shape = 4, stroke = 1.3
    ) +
    ggplot2::scale_color_manual(values = status_colors, name = "Recovery Status",
                                 drop = FALSE) +
    ggplot2::labs(
      x = "Parameter value", y = NULL,
      title = title %||% "Parameter Recovery",
      subtitle = "\u25C6 posterior mean   \u00D7 true value   bar = CI at classification width"
    ) +
    th$theme

  p
}