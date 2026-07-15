## ============================================================
## bsst_plot_status()
## Categorical recovery_status tile/mosaic, faceted by parameter
## ============================================================
#' Plot Recovery Status Across Two Design Variables (All Tracked Parameters)
#'
#' Tile plot of the categorical \code{recovery_status} (not a collapsed
#' scalar) across two design variables, faceted by parameter. Answers:
#' "across my design space, where does each tracked parameter individually
#' fall into each qualitative recovery bucket?" — deliberately not
#' collapsed to a single stress score, since averaging/maxing across
#' parameters (as \code{stress_score} does) can hide that, say, beta
#' recovers everywhere while sigma fails only when n is small.
#'
#' @param x Output list from \code{\link{bsst_stress}} or \code{\link{bsst_bo}}.
#' @param x_var,y_var Character; the two design variables for the axes.
#' @param parameters Character vector of tracked parameters to facet over.
#'   Required — must be supplied explicitly.
#' @param fix_others Named list of fixed values for other varying design
#'   variables. Numeric variables use nearest-match.
#' @param palette One of \code{"dark_research"}, \code{"github_dark"},
#'   \code{"scientific"}.
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#' @export
bsst_plot_status <- function(x, x_var, y_var, parameters,
                              fix_others = NULL,
                              palette = "dark_research",
                              title = NULL) {

  if (missing(parameters) || is.null(parameters) || length(parameters) == 0) {
    stop("`parameters` must be supplied explicitly (one or more tracked parameter names).")
  }
  if (!"raw_table" %in% names(x)) stop("x must be the output of bsst_stress() or bsst_bo().")

  df <- x$raw_table
  for (v in c(x_var, y_var)) if (!v %in% names(df)) stop("'", v, "' not found in x$raw_table.")

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

  df <- df[df$parameter %in% parameters, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows found for parameters: ", paste(parameters, collapse = ", "))

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

  if (any(duplicated(df[, c(x_var, y_var, "parameter")]))) {
    stop("Multiple rows share the same (", x_var, ", ", y_var, ", parameter) after filtering. ",
         "bsst_plot_status() expects a single gridded design (e.g. bsst_stress() full_factorial/lhs ",
         "output); it is not well-defined for scattered/duplicate designs like raw bsst_bo() output.")
  }

  df$recovery_status <- factor(df$recovery_status, levels = .bsst_status_levels)

  ggtheme <- .bsst_theme_ggplot(palette)
  pal <- ggtheme$palette
  status_colors <- setNames(pal$ordinal, .bsst_status_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = recovery_status)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ parameter) +
    ggplot2::scale_fill_manual(values = status_colors, name = "Recovery Status", drop = FALSE) +
    ggplot2::labs(
      x = x_var, y = y_var,
      title = title %||% "Recovery Status by Parameter",
      subtitle = if (length(other_varying) > 0) {
        paste0("Fixed: ", paste(sprintf("%s = %s", names(fix_others), unlist(fix_others)), collapse = ", "))
      } else NULL
    ) + ggtheme$theme

  p
}