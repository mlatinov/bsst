## ============================================================
## bsst_plot_surrogate()
## Fitted GP surface from bsst_bo() output — 2D heatmap or 3D surface
## ============================================================
#' Plot the Fitted GP Surrogate Surface from a BO Search
#'
#' Visualizes the Gaussian process surrogate fitted by \code{\link{bsst_bo}}
#' over the continuous search box, evaluated on a fine grid. Unlike
#' \code{bsst_plot_stress_2d()}/\code{_3d()} (which only show actually
#' observed points), this shows the surrogate's model-based interpolation
#' across the entire space, including regions never directly simulated —
#' answering "based on everything learned so far, what does the full error
#' surface look like?"
#'
#' @param x Output list from \code{\link{bsst_bo}}.
#' @param x_var,y_var Character; the two search-box dimensions to plot.
#'   Must both be present in \code{x$bounds}. If the search box had more
#'   than two dimensions, the remaining ones are fixed via \code{fix_others}
#'   (as fractions of their 0,1 normalized range, or in original units —
#'   see \code{fix_others}).
#' @param fix_others Named list of fixed values (in original design-variable
#'   units, not normalized) for any search dimensions besides
#'   \code{x_var}/\code{y_var}. Required if the search box had more than 2
#'   dimensions.
#' @param resolution Grid resolution per axis for surrogate evaluation.
#'   Default 50 (i.e. a 50x50 grid for 2D).
#' @param show_uncertainty Logical; if \code{TRUE}, also returns a second
#'   panel/plot showing the surrogate's predictive standard deviation
#'   instead of its mean, so regions of high uncertainty (under-explored)
#'   are visible alongside regions of high predicted stress. Only supported
#'   for \code{render = "surface_2d"}, returned as a list of two plots in
#'   that case.
#' @param render One of \code{"surface_2d"} (ggplot2 heatmap of the GP
#'   mean) or \code{"surface_3d"} (plotly 3D surface).
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object, a plotly object, or (if
#'   \code{show_uncertainty = TRUE} and \code{render = "surface_2d"}) a
#'   list with elements \code{mean} and \code{sd}.
#' @export
bsst_plot_surrogate <- function(x, x_var, y_var,
                                 fix_others = NULL,
                                 resolution = 50,
                                 show_uncertainty = FALSE,
                                 render = c("surface_2d", "surface_3d"),
                                 palette = "dark_research",
                                 title = NULL) {

  render <- match.arg(render)

  if (!"surrogate" %in% names(x) || is.null(x$surrogate)) {
    stop("x must be the output of bsst_bo() with a successfully fitted surrogate.")
  }
  if (missing(x_var) || missing(y_var)) stop("x_var and y_var must both be supplied.")

  bounds <- x$bounds
  all_dims <- names(bounds)
  if (!x_var %in% all_dims) stop("x_var '", x_var, "' not found in x$bounds.")
  if (!y_var %in% all_dims) stop("y_var '", y_var, "' not found in x$bounds.")

  other_dims <- setdiff(all_dims, c(x_var, y_var))
  if (length(other_dims) > 0) {
    missing_fix <- setdiff(other_dims, names(fix_others))
    if (length(missing_fix) > 0) {
      stop("The search box has other dimensions besides '", x_var, "' and '", y_var, "': ",
           paste(other_dims, collapse = ", "),
           ". You must supply `fix_others` (in original units) covering: ",
           paste(missing_fix, collapse = ", "))
    }
  }

  ## ---- build the evaluation grid in original units, then normalize ----
  x_seq <- seq(bounds[[x_var]][1], bounds[[x_var]][2], length.out = resolution)
  y_seq <- seq(bounds[[y_var]][1], bounds[[y_var]][2], length.out = resolution)
  grid_orig <- expand.grid(x_val = x_seq, y_val = y_seq)
  names(grid_orig) <- c(x_var, y_var)

  for (nm in other_dims) grid_orig[[nm]] <- fix_others[[nm]]

  grid_unit <- as.data.frame(matrix(nrow = nrow(grid_orig), ncol = length(all_dims)))
  names(grid_unit) <- all_dims
  for (nm in all_dims) {
    grid_unit[[nm]] <- (grid_orig[[nm]] - bounds[[nm]][1]) / (bounds[[nm]][2] - bounds[[nm]][1])
  }

  pred <- predict(x$surrogate, newdata = grid_unit, type = "UK", checkNames = FALSE)
  grid_orig$mean <- pred$mean
  grid_orig$sd   <- pred$sd

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette

  if (render == "surface_2d") {
    p_mean <- ggplot2::ggplot(grid_orig, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = mean)) +
      ggplot2::geom_raster(interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(colors = pal$sequential, name = "Predicted\nstress score") +
      ggplot2::labs(x = x_var, y = y_var,
                    title = title %||% "GP Surrogate \u2014 Predicted Stress Surface",
                    subtitle = if (length(other_dims) > 0) {
                      paste0("Fixed: ", paste(sprintf("%s = %s", names(fix_others), unlist(fix_others)), collapse = ", "))
                    } else NULL) +
      ggtheme$theme

    if (!show_uncertainty) return(p_mean)

    p_sd <- ggplot2::ggplot(grid_orig, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = sd)) +
      ggplot2::geom_raster(interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(colors = pal$sequential, name = "Predictive SD\n(uncertainty)") +
      ggplot2::labs(x = x_var, y = y_var, title = "GP Surrogate \u2014 Predictive Uncertainty") +
      ggtheme$theme

    return(list(mean = p_mean, sd = p_sd))
  }

  ## ---- surface_3d (plotly) ----
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for render = 'surface_3d'.")
  }
  layout_theme <- .bsst_theme_plotly_layout(palette)
  z_mat <- matrix(grid_orig$mean, nrow = length(x_seq), ncol = length(y_seq))
  colorscale_rgb <- grDevices::colorRampPalette(pal$sequential)(256)

  p <- plotly::plot_ly(
    x = ~y_seq, y = ~x_seq, z = ~z_mat, type = "surface",
    colors = colorscale_rgb, colorbar = list(title = "Predicted\nstress score")
  )
  p <- plotly::layout(
    p,
    title = title %||% "GP Surrogate \u2014 Predicted Stress Surface",
    scene = list(
      xaxis = c(list(title = x_var), layout_theme$xaxis),
      yaxis = c(list(title = y_var), layout_theme$yaxis),
      zaxis = list(title = "Predicted stress score"),
      bgcolor = layout_theme$plot_bgcolor
    ),
    paper_bgcolor = layout_theme$paper_bgcolor, font = layout_theme$font
  )
  p
}