## ============================================================
## bsst_plot_diagnostics()
## Compact pass/fail diagnostics panel
## ============================================================
#' Plot Convergence Diagnostics Summary
#'
#' Displays a compact pass/fail panel of sampler diagnostics (max R-hat,
#' min bulk/tail ESS, divergences) against the thresholds used in
#' \code{\link{bsst_recover}}, so a recovery plot is never read in
#' isolation from whether the underlying fit actually converged.
#'
#' @param x Output list from \code{\link{bsst_recover}}.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param thresholds Named list of thresholds to compare against; defaults
#'   to the thresholds stored in \code{x} if not supplied.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_diagnostics <- function(x, palette = "dark_research",
                                   thresholds = NULL, title = NULL) {

  if (!"diagnostics" %in% names(x)) {
    stop("x must be the output of bsst_recover().")
  }
  diag <- x$diagnostics
  th <- .bsst_theme_ggplot(palette)

  thresholds <- thresholds %||% list(
    rhat_max = 1.01, ess_bulk_min = 400, ess_tail_min = 400, divergences_max = 0
  )

  rows <- data.frame(
    metric = c("Max R-hat", "Min ESS (bulk)", "Min ESS (tail)", "Divergences"),
    value = c(diag$rhat_max, diag$ess_bulk_min, diag$ess_tail_min, diag$divergences),
    pass = c(
      diag$rhat_max <= thresholds$rhat_max,
      diag$ess_bulk_min >= thresholds$ess_bulk_min,
      diag$ess_tail_min >= thresholds$ess_tail_min,
      diag$divergences <= thresholds$divergences_max
    ),
    stringsAsFactors = FALSE
  )
  rows$metric <- factor(rows$metric, levels = rev(rows$metric))
  rows$label <- sprintf("%s: %.3f", rows$metric, rows$value)
  rows$status_color <- ifelse(rows$pass, th$palette$ordinal[1], th$palette$accent)

  p <- ggplot2::ggplot(rows, ggplot2::aes(x = 1, y = metric)) +
    ggplot2::geom_tile(ggplot2::aes(fill = pass), width = 0.95, height = 0.85, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label), color = th$palette$text_main,
                        hjust = 0.5, size = 3.6) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = th$palette$ordinal[1],
                                           `FALSE` = th$palette$accent)) +
    ggplot2::labs(
      title = title %||% "Convergence Diagnostics",
      subtitle = if (diag$converged) "All thresholds met" else "One or more thresholds violated",
      x = NULL, y = NULL
    ) +
    th$theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank()
    )

  p
}