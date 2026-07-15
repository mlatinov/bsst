## ============================================================
## bsst_plot_stress_2d()
## ============================================================
#' Plot Recovery/Stress Metric Across Two Design Variables (Heatmap)
#'
#' Visualizes a recovery metric as a 2D heatmap across two simulation
#' design variables, for a single tracked parameter. Answers: "is there an
#' interaction between two simulation conditions that jointly break
#' recovery — e.g. small sample size is fine unless measurement error is
#' also high?"
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var,y_var Character; names of the two design variables for the
#'   horizontal/vertical axes.
#' @param parameter Character; single tracked parameter name. Required.
#' @param fill_var One of \code{"std_error"}, \code{"error"},
#'   \code{"stress_score"}. Default \code{"stress_score"}.
#' @param fix_others Named list of fixed values for any other varying
#'   design variables, required whenever more than \code{x_var}/\code{y_var}
#'   vary in the underlying design. Numeric variables use nearest-match;
#'   non-numeric require an exact match.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param engine One of \code{"ggplot2"}, \code{"plotly"}.
#' @param show_nonconverged_marker Logical; if \code{TRUE} (default), tiles
#'   corresponding to a non-converged fit are marked with a small "x"
#'   overlay rather than left visually identical to converged tiles.
#' @param title Optional plot title.
#'
#' @return A ggplot2 or plotly object depending on \code{engine}.
#' @export
bsst_plot_stress_2d <- function(x, x_var, y_var, parameter,
                                 fill_var = c("stress_score", "std_error", "error"),
                                 fix_others = NULL,
                                 palette = "dark_research",
                                 engine = c("ggplot2", "plotly"),
                                 show_nonconverged_marker = TRUE,
                                 title = NULL) {

  fill_var <- match.arg(fill_var)
  engine <- match.arg(engine)

  if (missing(parameter) || is.null(parameter) || length(parameter) != 1) {
    stop("`parameter` must be supplied explicitly as a single parameter name.")
  }
  if (!"raw_table" %in% names(x)) {
    stop("x must be the output of bsst_stress() or bsst_bo().")
  }

  df <- x$raw_table
  for (v in c(x_var, y_var, fill_var)) {
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

  ## ---- for bsst_bo() output: points are scattered, not gridded ----
  ## detect duplicates at same (x_var, y_var) after filtering and aggregate defensively
  if (any(duplicated(df[, c(x_var, y_var)]))) {
    warning("Multiple rows share the same (", x_var, ", ", y_var, ") after filtering ",
            "(common with bsst_bo() output, whose points aren't gridded). Averaging ",
            "fill_var within duplicate cells.")
    df <- stats::aggregate(
      df[[fill_var]], by = list(x = df[[x_var]], y = df[[y_var]]), FUN = mean, na.rm = TRUE
    )
    names(df) <- c(x_var, y_var, fill_var)
    df$converged <- TRUE  # can't meaningfully preserve per-cell convergence after aggregation
  }

  pal <- .bsst_get_palette(palette)

  fill_label <- switch(fill_var,
    std_error = "Standardized error", error = "Error", stress_score = "Stress score")

  is_signed <- fill_var %in% c("std_error", "error")

  if (engine == "ggplot2") {
    ggtheme <- .bsst_theme_ggplot(palette)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
      ggplot2::geom_tile()

    if (is_signed) {
      lim <- max(abs(df[[fill_var]]), na.rm = TRUE)
      p <- p + ggplot2::scale_fill_gradientn(
        colors = pal$diverging, limits = c(-lim, lim), name = fill_label
      )
    } else {
      p <- p + ggplot2::scale_fill_gradientn(colors = pal$sequential, name = fill_label)
    }

    if (show_nonconverged_marker && "converged" %in% names(df)) {
      bad <- df[!df$converged, , drop = FALSE]
      if (nrow(bad) > 0) {
        p <- p + ggplot2::geom_point(data = bad, shape = 4, color = pal$text_main, size = 2)
      }
    }

    p <- p +
      ggplot2::labs(
        x = x_var, y = y_var,
        title = title %||% paste0(fill_label, " \u2014 ", parameter),
        subtitle = if (length(other_varying) > 0) {
          paste0("Fixed: ", paste(sprintf("%s = %s", names(fix_others), unlist(fix_others)), collapse = ", "))
        } else NULL
      ) +
      ggtheme$theme

    return(p)
  }

  ## ---- plotly engine ----
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for engine = 'plotly'. Install it or use engine = 'ggplot2'.")
  }

  layout_theme <- .bsst_theme_plotly_layout(palette)
  colorscale <- if (is_signed) {
    lim <- max(abs(df[[fill_var]]), na.rm = TRUE)
    list(colors = pal$diverging, zmin = -lim, zmax = lim)
  } else {
    list(colors = pal$sequential, zmin = min(df[[fill_var]], na.rm = TRUE),
         zmax = max(df[[fill_var]], na.rm = TRUE))
  }

  p <- plotly::plot_ly(
    data = df, x = ~get(x_var), y = ~get(y_var), z = ~get(fill_var),
    type = "heatmap",
    colors = grDevices::colorRampPalette(colorscale$colors)(256),
    zmin = colorscale$zmin, zmax = colorscale$zmax,
    hovertemplate = paste0(x_var, ": %{x}<br>", y_var, ": %{y}<br>", fill_label, ": %{z}<extra></extra>")
  )

  p <- plotly::layout(
    p,
    title = title %||% paste0(fill_label, " \u2014 ", parameter),
    xaxis = c(list(title = x_var), layout_theme$xaxis),
    yaxis = c(list(title = y_var), layout_theme$yaxis),
    paper_bgcolor = layout_theme$paper_bgcolor,
    plot_bgcolor  = layout_theme$plot_bgcolor,
    font = layout_theme$font
  )

  p
}