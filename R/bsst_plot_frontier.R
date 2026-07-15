## ============================================================
## bsst_plot_frontier()
## ============================================================
#' Plot the Acceptable/Unacceptable Recovery Frontier
#'
#' Draws the boundary in the design space where model recovery crosses
#' from acceptable to unacceptable, for a single tracked parameter.
#' Supports 1D (a single design variable, boundary shown as a reference
#' line/rug) and 2D (two design variables, boundary shown as a contour on
#' a heatmap). Answers the decision-relevant question: "what is the
#' minimum sample size (or maximum tolerable measurement error, etc.) for
#' acceptable inference?"
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var Character; design variable for the x-axis (1D or 2D mode).
#' @param y_var Character; design variable for the y-axis. If \code{NULL}
#'   (default), the function runs in 1D mode using \code{x_var} only.
#' @param parameter Character; single tracked parameter name. Required.
#' @param metric One of \code{"std_error"}, \code{"stress_score"}, used to
#'   evaluate a numeric \code{threshold} against. Ignored if
#'   \code{threshold} is a recovery-status string.
#' @param threshold Either a single number (interpreted against
#'   \code{metric}; points/regions with \code{abs(metric) > threshold} for
#'   \code{std_error}, or \code{metric > threshold} for \code{stress_score},
#'   are classified unacceptable) or a character string naming a
#'   \code{recovery_status} level (one of \code{"Recovered"},
#'   \code{"Recovered (wide)"}, \code{"Borderline"}, \code{"Not Recovered"},
#'   \code{"Severely Not Recovered"}) — in the latter case, the frontier is
#'   drawn at the boundary between that status and the next-worse one.
#' @param fix_others Named list of fixed values for any other varying
#'   design variables (1D mode: all besides \code{x_var}; 2D mode: all
#'   besides \code{x_var}/\code{y_var}). Numeric variables use
#'   nearest-match.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_frontier <- function(x, x_var, y_var = NULL, parameter,
                                metric = c("std_error", "stress_score"),
                                threshold,
                                fix_others = NULL,
                                palette = "dark_research",
                                title = NULL) {

  metric <- match.arg(metric)
  is_2d <- !is.null(y_var)

  if (missing(parameter) || is.null(parameter) || length(parameter) != 1) {
    stop("`parameter` must be supplied explicitly as a single parameter name.")
  }
  if (missing(threshold)) {
    stop("`threshold` must be supplied: either a number (evaluated against `metric`) ",
         "or a recovery_status level string (e.g. 'Not Recovered').")
  }
  if (!"raw_table" %in% names(x)) stop("x must be the output of bsst_stress() or bsst_bo().")

  df <- x$raw_table
  dvars <- .bsst_get_design_vars(x)
  if (!x_var %in% dvars$varying) stop("x_var '", x_var, "' is not a varying design variable.")
  if (is_2d && !y_var %in% dvars$varying) stop("y_var '", y_var, "' is not a varying design variable.")

  used_vars <- if (is_2d) c(x_var, y_var) else x_var
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

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette

  ## ---- determine pass/fail per row, and a numeric score for contouring ----
  if (is.character(threshold)) {
    if (!threshold %in% .bsst_status_levels) {
      stop("threshold string must be one of: ", paste(.bsst_status_levels, collapse = ", "))
    }
    status_rank <- setNames(seq_along(.bsst_status_levels), .bsst_status_levels)
    cutoff_rank <- status_rank[[threshold]]
    df$row_rank <- status_rank[df$recovery_status]
    df$acceptable <- df$row_rank <= cutoff_rank
    score_var <- "row_rank"
    score_label <- paste0("Recovery status \u2264 '", threshold, "'")
  } else if (is.numeric(threshold)) {
    if (!metric %in% names(df)) stop("metric '", metric, "' not found in x$raw_table.")
    df$eval_metric <- if (metric == "std_error") abs(df[[metric]]) else df[[metric]]
    df$acceptable <- df$eval_metric <= threshold
    score_var <- "eval_metric"
    score_label <- paste0(metric, " \u2264 ", threshold)
  } else {
    stop("threshold must be a number or a recovery_status string.")
  }

  df$status_label <- ifelse(df$acceptable, "Acceptable", "Unacceptable")

  if (!is_2d) {
    ## ---- 1D: line of the metric, shaded/marked region past the frontier ----
    df <- df[order(df[[x_var]]), ]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[score_var]])) +
      ggplot2::geom_line(color = pal$text_secondary, linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(color = status_label), size = 2.8) +
      ggplot2::scale_color_manual(values = c(Acceptable = pal$ordinal[1], Unacceptable = pal$accent),
                                   name = NULL)

    if (is.numeric(threshold)) {
      p <- p + ggplot2::geom_hline(yintercept = threshold, color = pal$text_main,
                                    linetype = "dashed", linewidth = 0.6)
    } else {
      p <- p + ggplot2::geom_hline(yintercept = cutoff_rank + 0.5, color = pal$text_main,
                                    linetype = "dashed", linewidth = 0.6)
    }

    p <- p + ggplot2::labs(
      x = x_var, y = score_label,
      title = title %||% paste0("Recovery Frontier \u2014 ", parameter),
      subtitle = paste0("Threshold: ", score_label)
    ) + ggtheme$theme

    return(p)
  }

  ## ---- 2D: heatmap of pass/fail with a contour line at the boundary ----
  if (any(duplicated(df[, c(x_var, y_var)]))) {
    warning("Multiple rows share the same (", x_var, ", ", y_var, ") after filtering; ",
            "averaging ", score_var, " within duplicate cells.")
    df <- stats::aggregate(df[[score_var]], by = list(xx = df[[x_var]], yy = df[[y_var]]), FUN = mean, na.rm = TRUE)
    names(df) <- c(x_var, y_var, score_var)
    df$acceptable <- if (is.numeric(threshold)) df[[score_var]] <= threshold else df[[score_var]] <= cutoff_rank
    df$status_label <- ifelse(df$acceptable, "Acceptable", "Unacceptable")
  }

  contour_break <- if (is.numeric(threshold)) threshold else cutoff_rank + 0.5

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_tile(ggplot2::aes(fill = status_label), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c(Acceptable = pal$ordinal[1], Unacceptable = pal$accent),
                                name = NULL) +
    ggplot2::geom_contour(ggplot2::aes(z = .data[[score_var]]), breaks = contour_break,
                           color = pal$text_main, linewidth = 0.9) +
    ggplot2::labs(
      x = x_var, y = y_var,
      title = title %||% paste0("Recovery Frontier \u2014 ", parameter),
      subtitle = paste0("Threshold: ", score_label)
    ) + ggtheme$theme

  p
}