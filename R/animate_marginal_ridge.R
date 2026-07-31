#' Animate a temporal posterior distribution in 3D (ODE, SDE, or any model)
#'
#' Renders posterior draws of a time-evolving quantity as an interactive 3D
#' scene: time on one axis, the quantity's value on another, and posterior
#' density on the third, colored on a fixed cold-to-hot thermal scale. Color
#' is normalized independently within each time slice, so the peak of every
#' time slice's own posterior reads as the hottest color and its own tails
#' read as coldest -- regardless of how that slice's absolute density
#' compares to any other slice's. Works identically for ODE closed-form fits,
#' SDE (OU/GBM/Euler-Maruyama) fits, or any model that produces posterior
#' draws over time -- the function only ever sees a long-format data frame of
#' (time, value, draw), so it is agnostic to which model family produced them.
#'
#' Two rendering modes:
#' \describe{
#'   \item{"surface"}{a smooth kernel-density ridge over time (the standard,
#'     ribbon-like view) -- easier to read at a glance.}
#'   \item{"cloud"}{every individual posterior draw plotted as a semi-transparent
#'     point, colored by its own local posterior density -- the raw,
#'     unsmoothed view. Slower to render for very large draw counts, so
#'     \code{n_draws_show} subsamples per time slice.}
#' }
#'
#' @param draws Long-format data.frame with required columns \code{time},
#'   \code{value}, \code{draw} (an integer draw index). Extract this from a
#'   Stan fit with e.g. \code{rstan::extract()} or \code{posterior::as_draws_df()}
#'   reshaped to long format -- one row per (draw, time) pair.
#' @param unit Optional value to filter \code{draws$unit} to a single unit
#'   (lake/subject/culture) if \code{draws} contains multiple. Required if
#'   the data frame has a \code{unit} column and more than one unit is present.
#' @param mode "surface" (smooth density ridge, default) or "cloud" (raw
#'   posterior draws as points, colored by local density).
#' @param theme One of "midnight" (dark, cool), "sunrise" (warm, light), or
#'   "mono" (grayscale, minimal). Purely cosmetic -- pick whichever reads best
#'   for the audience/output medium.
#' @param observed Optional data.frame with columns \code{time}, \code{value}
#'   -- e.g. the actual simulated/observed data -- overlaid as a distinct
#'   marker series on top of the posterior.
#' @param n_draws_show Max number of posterior draws rendered per time slice
#'   in "cloud" mode (subsampled for rendering performance). Ignored in
#'   "surface" mode.
#' @param bandwidth Kernel bandwidth for the density estimate. Defaults to
#'   1/20th of the value range if not supplied -- tune if the ridge/cloud
#'   coloring looks too spiky or too smoothed out.
#' @param value_range Optional c(min, max) for the value axis. Defaults to
#'   the observed range of \code{draws$value} with 5% padding.
#' @param title Optional plot title.
#'
#' @return A plotly htmlwidget. Assign it and print/View it, or pipe into
#'   \code{htmlwidgets::saveWidget()} to export a standalone HTML file.
#'
#' @examples
#' \dontrun{
#' # from a Stan fit's generated quantities (e.g. y_rep), reshaped long:
#' post <- as.data.frame(rstan::extract(fit, "y_rep")$y_rep)
#' names(post) <- paste0("t", seq_len(ncol(post)))
#' post$draw <- seq_len(nrow(post))
#' long <- tidyr::pivot_longer(post, -draw, names_to = "time",
#'                              names_prefix = "t", values_to = "value")
#' long$time <- as.numeric(long$time)
#'
#' animate_temporal_posterior(long, mode = "cloud", theme = "midnight")
#' animate_temporal_posterior(long, mode = "surface", theme = "sunrise",
#'                             observed = sim$data)
#' }
#'
#' @export
animate_temporal_posterior <- function(draws,
                                        unit = NULL,
                                        mode = c("surface", "cloud"),
                                        theme = c("midnight", "sunrise", "mono"),
                                        observed = NULL,
                                        n_draws_show = 400,
                                        bandwidth = NULL,
                                        value_range = NULL,
                                        title = NULL) {

  mode  <- match.arg(mode)
  theme <- match.arg(theme)

  stopifnot(all(c("time", "value", "draw") %in% names(draws)))

  if (!is.null(unit)) {
    stopifnot("unit" %in% names(draws))
    draws <- draws[draws$unit == unit, ]
  } else if ("unit" %in% names(draws) && length(unique(draws$unit)) > 1) {
    stop("`draws` contains multiple units -- pass `unit = <id>` to select one.")
  }

  times <- sort(unique(draws$time))

  if (is.null(value_range)) {
    value_range <- range(draws$value, na.rm = TRUE)
    pad <- diff(value_range) * 0.05
    value_range <- value_range + c(-pad, pad)
  }
  if (is.null(bandwidth)) {
    bandwidth <- diff(value_range) / 20
  }

  # one kernel density estimate per time slice -- reused for both modes:
  # "surface" turns these into the ridge; "cloud" uses them to color points
  # by their own local density.
  dens_by_time <- lapply(times, function(t) {
    vals <- draws$value[draws$time == t]
    stats::density(vals, bw = bandwidth,
                    from = value_range[1], to = value_range[2], n = 256)
  })
  names(dens_by_time) <- as.character(times)

  # theme controls page chrome only (background / ink / observed-point accent) --
  # density coloring below is a fixed perceptual scale, independent of theme,
  # so "how dense is this region" always reads the same way regardless of
  # which page theme is active.
  palette <- switch(theme,
    midnight = list(bg = "#0b0c10", ink = "#c9d6e3", accent = "#ffd166"),
    sunrise  = list(bg = "#fff8f0", ink = "#4a3527", accent = "#2e6f95"),
    mono     = list(bg = "#f4f4f2", ink = "#333333", accent = "#2e6f95")
  )

  # cold (low density) -> vivid hot (peak density): navy -> violet -> magenta
  # -> orange -> red. Fixed across themes so density always reads as "heat".
  colorscale <- list(
    list(0,    "#0d1b4c"),
    list(0.35, "#5b2a86"),
    list(0.6,  "#b5308f"),
    list(0.8,  "#e8622f"),
    list(1,    "#ff2b2b")
  )

  # per-time-slice normalization: each time's own peak maps to 1 (hottest
  # color), each time's own lowest density maps to 0 (coldest) -- so "hot"
  # always means "the peak of THIS time slice's posterior", not a value
  # compared across slices. Absolute density magnitude differs slice to
  # slice for reasons unrelated to the story we want to tell (bandwidth,
  # local sample density), so we deliberately don't compare across time.
  dens_max_by_time <- vapply(dens_by_time, function(d) max(d$y), numeric(1))

  if (mode == "cloud") {

    draws$density <- mapply(function(t, v) {
      d <- dens_by_time[[as.character(t)]]
      stats::approx(d$x, d$y, xout = v)$y
    }, draws$time, draws$value)

    draws$density_norm <- mapply(function(t, dens) {
      dens / dens_max_by_time[as.character(t)]
    }, draws$time, draws$density)

    draws_show <- do.call(rbind, lapply(split(draws, draws$time), function(df) {
      if (nrow(df) > n_draws_show) df <- df[sample(nrow(df), n_draws_show), ]
      df
    }))

    p <- plotly::plot_ly(
      data = draws_show, x = ~time, y = ~value, z = ~density,
      type = "scatter3d", mode = "markers",
      marker = list(size = 3, opacity = 0.35,
                    color = ~density_norm, colorscale = colorscale,
                    cmin = 0, cmax = 1, showscale = FALSE)
    )

  } else {

    z_mat <- sapply(dens_by_time, function(d) d$y)   # value-grid (rows) x time (cols)
    color_mat <- sweep(z_mat, 2, dens_max_by_time, `/`)  # each column normalized by its own peak
    p <- plotly::plot_ly(
      x = times, y = dens_by_time[[1]]$x, z = z_mat,
      surfacecolor = color_mat,
      type = "surface", colorscale = colorscale,
      cmin = 0, cmax = 1,
      showscale = FALSE, opacity = 0.96,
      contours = list(z = list(show = FALSE))
    )

  }

  if (!is.null(observed)) {
    stopifnot(all(c("time", "value") %in% names(observed)))
    obs_dens <- mapply(function(t, v) {
      d <- dens_by_time[[as.character(t)]]
      if (is.null(d)) return(0)
      stats::approx(d$x, d$y, xout = v)$y
    }, observed$time, observed$value)

    p <- plotly::add_trace(
      p, data = observed, x = ~time, y = ~value, z = obs_dens,
      type = "scatter3d", mode = "markers",
      marker = list(size = 4, color = palette$accent, symbol = "diamond"),
      inherit = FALSE, name = "observed"
    )
  }

  p <- plotly::layout(p,
    paper_bgcolor = palette$bg, plot_bgcolor = palette$bg,
    font = list(color = palette$ink),
    title = title,
    scene = list(
      xaxis = list(title = "time", color = palette$ink, backgroundcolor = palette$bg,
                   gridcolor = palette$ink),
      yaxis = list(title = "value", color = palette$ink, backgroundcolor = palette$bg,
                   gridcolor = palette$ink),
      zaxis = list(title = "posterior density", color = palette$ink,
                   backgroundcolor = palette$bg, gridcolor = palette$ink),
      bgcolor = palette$bg,
      camera = list(eye = list(x = 1.7, y = 1.7, z = 0.85))
    )
  )

  p
}