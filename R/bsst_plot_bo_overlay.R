## ============================================================
## bsst_plot_bo_overlay()
## Evaluated points colored by origin, over the surrogate surface
## ============================================================
#' Plot Evaluated Points Over the Surrogate Surface, Colored by Origin
#'
#' Overlays the actually-evaluated design points from \code{\link{bsst_bo}}
#' (colored by whether they came from a warm start, the random initial
#' design, or adaptive BO batches) on top of the fitted surrogate surface.
#' Answers: "did the optimizer actually concentrate its search around the
#' worst region, or is it still scattered / exploring broadly?"
#'
#' @param x Output list from \code{\link{bsst_bo}}.
#' @param x_var,y_var Character; the two search-box dimensions to plot.
#' @param fix_others Named list of fixed values for other search
#'   dimensions, as in \code{\link{bsst_plot_surrogate}}.
#' @param resolution Grid resolution for the underlying surrogate surface.
#' @param origin Optional character vector, same length as the number of
#'   evaluated points, labeling each point's origin (e.g.
#'   \code{c("warm_start", "warm_start", "init", "bo", ...)}). If
#'   \code{NULL} (default), all points are labeled \code{"evaluated"} since
#'   \code{bsst_bo()}'s current return object does not track per-point
#'   origin — see note below.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_bo_overlay <- function(x, x_var, y_var,
                                  fix_others = NULL,
                                  resolution = 50,
                                  origin = NULL,
                                  palette = "dark_research",
                                  title = NULL) {

  base_plot <- bsst_plot_surrogate(
    x, x_var = x_var, y_var = y_var, fix_others = fix_others,
    resolution = resolution, show_uncertainty = FALSE,
    render = "surface_2d", palette = palette, title = title %||% "BO Search Overlay"
  )

  df <- x$raw_table
  pts <- df[!duplicated(df$design_point_id), c("design_point_id", x_var, y_var)]

  if (is.null(origin)) {
    message("`origin` not supplied — bsst_bo()'s current return object does not track per-point ",
            "origin (warm start vs. init vs. BO batch), so all points are shown as a single group. ",
            "To distinguish origins, either supply `origin` manually (matching evaluation order in ",
            "x$raw_table's unique design_point_id order) or extend bsst_bo() to record and return it.")
    pts$origin <- "evaluated"
  } else {
    if (length(origin) != nrow(pts)) {
      stop("`origin` length (", length(origin), ") must match the number of unique evaluated points (",
           nrow(pts), ").")
    }
    pts$origin <- origin
  }

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette
  origin_levels <- unique(pts$origin)
  origin_colors <- setNames(
    grDevices::colorRampPalette(c(pal$text_main, pal$accent))(length(origin_levels)),
    origin_levels
  )

  p <- base_plot +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], color = origin),
                         size = 2.2, inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = origin_colors, name = "Origin")

  p
}