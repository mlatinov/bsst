## ============================================================
## bsst_plot_bo_trace()
## Running best-found stress score vs. BO iteration
## ============================================================
#' Plot Bayesian Optimization Search Progress
#'
#' Plots the running-best (cumulative max) stress score found by
#' \code{\link{bsst_bo}} against evaluation order, with individual batch
#' proposals shown as fainter points. Answers: "has the search converged
#' on a worst case, or is it still finding progressively worse regions —
#' do I need more iterations?"
#'
#' @param x Output list from \code{\link{bsst_bo}}.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_bo_trace <- function(x, palette = "dark_research", title = NULL) {

  if (!"raw_table" %in% names(x)) stop("x must be the output of bsst_bo().")

  df <- x$raw_table
  if (!"design_point_id" %in% names(df) || !"stress_score" %in% names(df)) {
    stop("x$raw_table must contain design_point_id and stress_score columns.")
  }

  pts <- df[!duplicated(df$design_point_id), c("design_point_id", "stress_score")]
  pts <- pts[order(pts$design_point_id), ]
  pts$evaluation_order <- seq_len(nrow(pts))
  pts$running_best <- cummax(pts$stress_score)

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = evaluation_order)) +
    ggplot2::geom_point(ggplot2::aes(y = stress_score), color = pal$text_secondary,
                         alpha = 0.5, size = 1.8) +
    ggplot2::geom_step(ggplot2::aes(y = running_best), color = pal$accent, linewidth = 0.9) +
    ggplot2::labs(
      x = "Evaluation order", y = "Stress score",
      title = title %||% "BO Search Progress",
      subtitle = "Faint points = individual evaluations; line = running best found"
    ) + ggtheme$theme

  p
}