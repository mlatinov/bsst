#' Animate the joint posterior probability cloud in 3D across time (or any index)
#'
#' Renders the JOINT posterior geometry of a multivariate state
#' x(t) = (x1(t), x2(t), x3(t)) at each time point (or any ordering variable
#' -- dose, sample size, culture id, a sweep parameter). Each point is
#' colored by a *local multivariate* density estimate, rank-normalized
#' before mapping to color so the color scale isn't crushed by the heavy
#' right-skew typical of kNN density estimates (most points near zero, a
#' few outliers very high) -- without this, almost everything renders near
#' the dark end of the scale, which looks like "everything is black."
#'
#' Colors are precomputed as literal rgba strings in R (not handed to
#' plotly as a numeric array + colorscale) because numeric marker.color +
#' colorscale on scatter3d combined with animation frames is unreliable in
#' Plotly.js -- baking the color avoids that failure mode entirely.
#'
#' New draws are added each frame and OLD draws remain visible, fading
#' toward a floor opacity rather than disappearing -- a growing, trailing
#' cloud instead of a preloaded backdrop or a hard cut between frames.
#'
#' Honest limitation: Plotly.js's frame-to-frame tweening is much better
#' supported for 2D scatter than scatter3d. The trail/fade mechanic here
#' makes the lack of true 3D interpolation much less noticeable, but it
#' does not manufacture smooth motion between frames -- for genuinely
#' continuous-looking deformation, a volumetric isosurface built from a
#' 3D KDE grid (rendered as translucent nested density shells instead of
#' discrete points) is the more fundamental fix and a good next step.
#'
#' @param draws Long-format data.frame: one row per (draw, time), with
#'   columns for the ordering variable, a draw id, and each state dimension.
#' @param dims Character vector of 2 or 3 column names in `draws` giving
#'   the state-space coordinates (e.g. c("x1","x2","x3")).
#' @param time_col Column name for the ordering/animation variable
#'   (default "time"). Can be any continuous sweep variable, not just time.
#' @param draw_col Column name for the draw index (default "draw").
#' @param density_method "knn" (fast, default) or "kde" (slower, more
#'   principled for small N via ks::kde).
#' @param k Neighbors used for "knn" density (default 15).
#' @param color_transform "rank" (default -- percentile rank, robust to
#'   skew), "log", or "identity". "rank" is almost always the right choice;
#'   use "identity" only if you specifically want raw density magnitude
#'   preserved in the color mapping.
#' @param n_draws_show Number of distinct draw ids used, same subset at
#'   every time point (object constancy).
#' @param theme "midnight", "sunrise", or "mono".
#' @param observed Optional data.frame with the same dims + time_col
#'   columns (e.g. the true/simulated trajectory), accumulated the same way
#'   as the posterior cloud but at full opacity (ground truth doesn't fade).
#' @param show_ridge Logical; show the mean-coordinate trajectory across
#'   all time points as a persistent line, always fully drawn.
#' @param trail_length Number of past time points to keep visible at once.
#'   `Inf` (default) keeps everything ever shown (matches "leave the
#'   original points there"). For long simulations this can produce a lot
#'   of rows (roughly T^2/2 * n_draws_show in the worst case) -- set a
#'   finite value (e.g. 10-15) if the browser gets sluggish.
#' @param trail_decay Per-frame opacity decay for older points (default
#'   0.8). Lower = fades faster.
#' @param trail_floor Minimum opacity older points fade to (default 0.15) --
#'   never fully invisible.
#' @param current_alpha Opacity of the newest frame's points (default 0.95).
#' @param frame_ms Milliseconds each frame is displayed during Play.
#' @param transition_ms Milliseconds for the tween between frames (mainly
#'   effective in 2D; scatter3d transition support is limited by Plotly.js).
#'
#' @return A plotly htmlwidget with a working Play/Pause button and a
#'   slider stepping through `time_col`. Assign and print, or
#'   htmlwidgets::saveWidget() to export a standalone HTML file.
#'
#' @examples
#' \dontrun{
#' animate_posterior_cloud(draws, dims = c("x1","x2","x3"),
#'                          theme = "midnight", trail_length = 12)
#' }
#'
#' @export
animate_posterior_cloud <- function(draws,
                                     dims = c("x1", "x2", "x3"),
                                     time_col = "time",
                                     draw_col = "draw",
                                     density_method = c("knn", "kde"),
                                     k = 15,
                                     color_transform = c("rank", "log", "identity"),
                                     n_draws_show = 1200,
                                     theme = c("midnight", "sunrise", "mono"),
                                     observed = NULL,
                                     show_ridge = TRUE,
                                     trail_length = Inf,
                                     trail_decay = 0.8,
                                     trail_floor = 0.15,
                                     current_alpha = 0.95,
                                     frame_ms = 500,
                                     transition_ms = 300) {

  density_method  <- match.arg(density_method)
  color_transform <- match.arg(color_transform)
  theme <- match.arg(theme)
  stopifnot(length(dims) %in% c(2, 3))
  stopifnot(all(c(dims, time_col, draw_col) %in% names(draws)))

  is3d  <- length(dims) == 3
  times <- sort(unique(draws[[time_col]]))
  T_n   <- length(times)

  palette <- switch(theme,
    midnight = list(bg = "#000004", ink = "#e8e2d5", ridge = "#4fd1c5", accent = "#4fd1c5"),
    sunrise  = list(bg = "#fff8f0", ink = "#4a3527", ridge = "#2e6f95", accent = "#2e6f95"),
    mono     = list(bg = "#0d0d0d", ink = "#dddddd", ridge = "#ffffff", accent = "#ffffff")
  )

  # magma-style ramp: stays legible on black because the low end is a
  # visible purple, not near-black navy
  magma_ramp <- grDevices::colorRampPalette(
    c("#000004", "#51127c", "#b73779", "#fc8961", "#fcfdbf")
  )(256)
  density_to_hex <- function(x01) {
    idx <- pmax(1, pmin(256, round(x01 * 255) + 1))
    magma_ramp[idx]
  }
  hex_to_rgba <- function(hex, alpha) {
    rgb <- grDevices::col2rgb(hex)
    sprintf("rgba(%d,%d,%d,%.3f)", rgb[1, ], rgb[2, ], rgb[3, ], alpha)
  }
  transform_density <- function(dens, method) {
    if (method == "rank") {
      rank(dens, ties.method = "average") / length(dens)
    } else if (method == "log") {
      ld <- log(dens - min(dens) + 1e-8)
      (ld - min(ld)) / (diff(range(ld)) + 1e-8)
    } else {
      (dens - min(dens)) / (diff(range(dens)) + 1e-8)
    }
  }

  all_ids  <- unique(draws[[draw_col]])
  keep_ids <- if (length(all_ids) > n_draws_show) sample(all_ids, n_draws_show) else all_ids
  draws <- draws[draws[[draw_col]] %in% keep_ids, , drop = FALSE]

  knn_density <- function(mat, k) {
    n <- nrow(mat); k <- min(k, n - 1)
    nn <- FNN::get.knn(mat, k = k)
    1 / (rowMeans(nn$nn.dist) ^ ncol(mat) + 1e-8)
  }
  kde_density <- function(mat) {
    fit <- ks::kde(mat)
    ks::predict.kde(fit, x = mat)
  }

  frames_data <- lapply(times, function(t) {
    df  <- draws[draws[[time_col]] == t, , drop = FALSE]
    mat <- as.matrix(df[, dims])
    df$density <- if (density_method == "knn") knn_density(mat, k) else kde_density(mat)
    df
  })

  pooled_density <- unlist(lapply(frames_data, `[[`, "density"))
  color01_all <- transform_density(pooled_density, color_transform)
  offset <- 0
  frames_data <- lapply(frames_data, function(df) {
    n <- nrow(df)
    df$color01 <- color01_all[(offset + 1):(offset + n)]
    df$hex <- density_to_hex(df$color01)
    offset <<- offset + n
    df
  })

  all_frames_df <- do.call(rbind, frames_data)

  axis_ranges <- lapply(dims, function(d) {
    v <- draws[[d]]; pad <- diff(range(v)) * 0.08
    range(v) + c(-pad, pad)
  })
  names(axis_ranges) <- dims

  ridge_df <- do.call(rbind, lapply(frames_data, function(df) {
    row <- as.data.frame(t(colMeans(df[, dims, drop = FALSE])))
    row[[time_col]] <- unique(df[[time_col]])
    row
  }))

  build_reveal <- function() {
    out <- vector("list", T_n)
    for (i in seq_len(T_n)) {
      j_start <- if (is.finite(trail_length)) max(1, i - trail_length + 1) else 1
      chunk <- do.call(rbind, lapply(j_start:i, function(j) {
        df <- frames_data[[j]]
        age <- i - j
        alpha <- if (age == 0) current_alpha else max(trail_floor, current_alpha * trail_decay ^ age)
        df$rgba <- hex_to_rgba(df$hex, alpha)
        df$point_id <- paste(df[[draw_col]], df[[time_col]], sep = "_")
        df$reveal_frame <- times[i]
        df
      }))
      out[[i]] <- chunk
    }
    do.call(rbind, out)
  }
  reveal_df <- build_reveal()

  # standardize dimension names to allow hardcoded ~dim1/~dim2/~dim3 formulas
  # below -- plotly's `frame` and `ids` arguments MUST be formulas bound to
  # a `data=` argument to be registered as animation drivers (this is what
  # builds the slider); passing raw vectors silently skips that
  # registration, which is why animation_slider() previously errored with
  # "attempt to select less than one element" -- there was no slider to find.
  reveal_df$dim1 <- reveal_df[[dims[1]]]
  reveal_df$dim2 <- reveal_df[[dims[2]]]
  if (is3d) reveal_df$dim3 <- reveal_df[[dims[3]]]

  ridge_df$dim1 <- ridge_df[[dims[1]]]
  ridge_df$dim2 <- ridge_df[[dims[2]]]
  if (is3d) ridge_df$dim3 <- ridge_df[[dims[3]]]

  p <- if (is3d) {
    plotly::plot_ly(data = reveal_df, x = ~dim1, y = ~dim2, z = ~dim3,
                     frame = ~reveal_frame, ids = ~point_id,
                     type = "scatter3d", mode = "markers",
                     marker = list(size = 3.5, color = reveal_df$rgba),
                     name = "posterior draws")
  } else {
    plotly::plot_ly(data = reveal_df, x = ~dim1, y = ~dim2,
                     frame = ~reveal_frame, ids = ~point_id,
                     type = "scatter", mode = "markers",
                     marker = list(size = 7, color = reveal_df$rgba),
                     name = "posterior draws")
  }

  # ridge is static (no frame column in ridge_df) -- inherit = FALSE stops
  # it from trying to inherit the base trace's frame mapping, which is
  # what caused the slider-registration error before
  if (show_ridge) {
    if (is3d) {
      p <- plotly::add_trace(p, data = ridge_df, x = ~dim1, y = ~dim2, z = ~dim3,
                              type = "scatter3d", mode = "lines+markers",
                              line = list(color = palette$ridge, width = 3),
                              marker = list(size = 3, color = palette$ridge),
                              name = "mean trajectory", showlegend = FALSE,
                              inherit = FALSE)
    } else {
      p <- plotly::add_trace(p, data = ridge_df, x = ~dim1, y = ~dim2,
                              type = "scatter", mode = "lines+markers",
                              line = list(color = palette$ridge, width = 3),
                              marker = list(size = 4, color = palette$ridge),
                              name = "mean trajectory", showlegend = FALSE,
                              inherit = FALSE)
    }
  }

  if (!is.null(observed)) {
    stopifnot(all(c(dims, time_col) %in% names(observed)))
    obs_reveal <- do.call(rbind, lapply(seq_along(times), function(i) {
      sub <- observed[observed[[time_col]] <= times[i], , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
      sub$reveal_frame <- times[i]
      sub
    }))
    obs_reveal$dim1 <- obs_reveal[[dims[1]]]
    obs_reveal$dim2 <- obs_reveal[[dims[2]]]
    if (is3d) obs_reveal$dim3 <- obs_reveal[[dims[3]]]

    if (is3d) {
      p <- plotly::add_trace(p, data = obs_reveal, x = ~dim1, y = ~dim2, z = ~dim3,
                              frame = ~reveal_frame,
                              type = "scatter3d", mode = "markers",
                              marker = list(size = 6, color = palette$accent, symbol = "diamond"),
                              name = "observed", inherit = FALSE)
    } else {
      p <- plotly::add_trace(p, data = obs_reveal, x = ~dim1, y = ~dim2,
                              frame = ~reveal_frame,
                              type = "scatter", mode = "markers",
                              marker = list(size = 10, color = palette$accent, symbol = "diamond"),
                              name = "observed", inherit = FALSE)
    }
  }

  p <- plotly::animation_opts(p, frame = frame_ms, transition = transition_ms,
                               easing = "cubic-in-out", redraw = TRUE)

  p <- plotly::layout(p, paper_bgcolor = palette$bg, plot_bgcolor = palette$bg,
                       font = list(color = palette$ink),
                       sliders = list(list(
                         currentvalue = list(prefix = paste0(time_col, ": "),
                                              font = list(color = palette$ink))
                       )))

  if (is3d) {
    p <- plotly::layout(p, scene = list(
      xaxis = list(title = dims[1], range = axis_ranges[[dims[1]]], color = palette$ink,
                   backgroundcolor = palette$bg, gridcolor = palette$ink),
      yaxis = list(title = dims[2], range = axis_ranges[[dims[2]]], color = palette$ink,
                   backgroundcolor = palette$bg, gridcolor = palette$ink),
      zaxis = list(title = dims[3], range = axis_ranges[[dims[3]]], color = palette$ink,
                   backgroundcolor = palette$bg, gridcolor = palette$ink),
      bgcolor = palette$bg
    ))
  } else {
    p <- plotly::layout(p,
      xaxis = list(title = dims[1], range = axis_ranges[[dims[1]]], color = palette$ink, gridcolor = palette$ink),
      yaxis = list(title = dims[2], range = axis_ranges[[dims[2]]], color = palette$ink, gridcolor = palette$ink)
    )
  }

  p
}