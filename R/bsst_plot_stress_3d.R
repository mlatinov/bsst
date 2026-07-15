## ============================================================
## bsst_plot_stress_3d()
## ============================================================
#' Plot Recovery/Stress Metric as a 3D Surface or Scatter (plotly)
#'
#' Visualizes a recovery metric across two design variables as a 3D
#' surface/scatter (fill/z-height = the metric), for a single tracked
#' parameter. Same underlying question as \code{\link{bsst_plot_stress_2d}}
#' — "is there an interaction between two simulation conditions that jointly
#' break recovery?" — but rendered with height as a third visual channel,
#' which some readers find easier to interpret for magnitude than color
#' alone. plotly-only, since ggplot2 has no native 3D support.
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var,y_var Character; the two design variables for the
#'   horizontal plane axes.
#' @param parameter Character; single tracked parameter name. Required.
#' @param z_var One of \code{"std_error"}, \code{"error"},
#'   \code{"stress_score"}. Default \code{"stress_score"}.
#' @param fix_others Named list of fixed values for any other varying
#'   design variables. Numeric variables use nearest-match.
#' @param render One of \code{"auto"}, \code{"surface"}, \code{"scatter"}.
#'   \code{"surface"} requires the data to form (or be interpolated onto) a
#'   regular grid — appropriate for \code{bsst_stress()} full-factorial
#'   output. \code{"scatter"} plots each evaluated point directly as a 3D
#'   point cloud with no interpolation — always valid, and the only sane
#'   choice for irregular designs like LHS or \code{bsst_bo()} output.
#'   \code{"auto"} (default) picks \code{"surface"} if the filtered data
#'   forms a regular grid, otherwise falls back to \code{"scatter"} with a
#'   message.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A plotly object.
#' @export
bsst_plot_stress_3d <- function(x, x_var, y_var, parameter,
                                 z_var = c("stress_score", "std_error", "error"),
                                 fix_others = NULL,
                                 render = c("auto", "surface", "scatter"),
                                 palette = "dark_research",
                                 title = NULL) {

  z_var <- match.arg(z_var)
  render <- match.arg(render)

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for bsst_plot_stress_3d(). Install it to use this function.")
  }
  if (missing(parameter) || is.null(parameter) || length(parameter) != 1) {
    stop("`parameter` must be supplied explicitly as a single parameter name.")
  }
  if (!"raw_table" %in% names(x)) {
    stop("x must be the output of bsst_stress() or bsst_bo().")
  }

  df <- x$raw_table
  for (v in c(x_var, y_var, z_var)) {
    if (!v %in% names(df)) stop("'", v, "' not found in x$raw_table.")
  }

  dvars <- .bsst_get_design_vars(x)
  if (!x_var %in% dvars$varying) stop("x_var '", x_var, "' is not a varying design variable.")
  if (!y_var %in% dvars$varying) stop("y_var '", y_var, "' is not a varying design variable.")

  other_varying <- setdiff(dvars$varying, c(x_var, y_var))
  if (length(other_varying) > 0) {
    missing_fix <- setdiff(other_varying, names(fix_others))
    if (length(missing_fix) > 0) {
      stop("The design has other varying variables besides '", x_var, "' and '", y_var, "': ",
           paste(other_varying, collapse = ", "),
           ". You must supply `fix_others` covering: ", paste(missing_fix, collapse = ", "))
    }
  }

  ## ---- filter to requested parameter ----
  df <- df[df$parameter == parameter, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows found for parameter '", parameter, "'.")

  ## ---- filter other varying dims (nearest-match for numeric) ----
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

  ## de-dupe identical (x_var, y_var) pairs by averaging z_var (same rationale as _2d)
  if (any(duplicated(df[, c(x_var, y_var)]))) {
    df <- stats::aggregate(df[[z_var]], by = list(xx = df[[x_var]], yy = df[[y_var]]), FUN = mean, na.rm = TRUE)
    names(df) <- c(x_var, y_var, z_var)
  }

  pal <- .bsst_get_palette(palette)
  z_label <- switch(z_var, std_error = "Standardized error", error = "Error", stress_score = "Stress score")
  is_signed <- z_var %in% c("std_error", "error")
  layout_theme <- .bsst_theme_plotly_layout(palette)

  ## ---- decide whether data forms a regular grid ----
  x_vals <- sort(unique(df[[x_var]]))
  y_vals <- sort(unique(df[[y_var]]))
  is_regular_grid <- nrow(df) == length(x_vals) * length(y_vals)

  if (render == "auto") {
    render <- if (is_regular_grid) "surface" else "scatter"
    if (render == "scatter") {
      message("Data does not form a regular grid over (", x_var, ", ", y_var,
              ") — using render = 'scatter'. This is expected for LHS or bsst_bo() output. ",
              "For an interpolated continuous surface from bsst_bo() output, use bsst_plot_surrogate() instead.")
    }
  } else if (render == "surface" && !is_regular_grid) {
    stop("render = 'surface' requires data forming a regular grid over (", x_var, ", ", y_var,
         "), but the filtered data does not (", nrow(df), " points for ", length(x_vals),
         " x ", length(y_vals), " unique values). Use render = 'scatter' or render = 'auto', ",
         "or use bsst_plot_surrogate() for an interpolated view of bsst_bo() output.")
  }

  color_scale <- if (is_signed) {
    lim <- max(abs(df[[z_var]]), na.rm = TRUE)
    list(colors = pal$diverging, zmin = -lim, zmax = lim)
  } else {
    list(colors = pal$sequential, zmin = min(df[[z_var]], na.rm = TRUE), zmax = max(df[[z_var]], na.rm = TRUE))
  }
  colorscale_rgb <- grDevices::colorRampPalette(color_scale$colors)(256)

  if (render == "surface") {
    z_mat <- matrix(NA_real_, nrow = length(x_vals), ncol = length(y_vals))
    for (i in seq_len(nrow(df))) {
      xi <- match(df[[x_var]][i], x_vals)
      yi <- match(df[[y_var]][i], y_vals)
      z_mat[xi, yi] <- df[[z_var]][i]
    }

    p <- plotly::plot_ly(
      x = ~y_vals, y = ~x_vals, z = ~z_mat, type = "surface",
      colors = colorscale_rgb, cmin = color_scale$zmin, cmax = color_scale$zmax,
      colorbar = list(title = z_label)
    )
  } else {
    p <- plotly::plot_ly(
      data = df, x = ~get(x_var), y = ~get(y_var), z = ~get(z_var),
      type = "scatter3d", mode = "markers",
      marker = list(
        size = 4,
        color = ~get(z_var),
        colorscale = list(seq(0, 1, length.out = length(colorscale_rgb)), colorscale_rgb),
        cmin = color_scale$zmin, cmax = color_scale$zmax,
        colorbar = list(title = z_label)
      ),
      hovertemplate = paste0(x_var, ": %{x}<br>", y_var, ": %{y}<br>", z_label, ": %{z}<extra></extra>")
    )
  }

  p <- plotly::layout(
    p,
    title = title %||% paste0(z_label, " \u2014 ", parameter),
    scene = list(
      xaxis = c(list(title = x_var), layout_theme$xaxis),
      yaxis = c(list(title = y_var), layout_theme$yaxis),
      zaxis = c(list(title = z_label), layout_theme$xaxis),
      bgcolor = layout_theme$plot_bgcolor
    ),
    paper_bgcolor = layout_theme$paper_bgcolor,
    font = layout_theme$font
  )

  p
}