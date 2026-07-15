library(ggplot2)

## ============================================================
## Palette registry
## ============================================================
.bsst_palettes <- list(
  dark_research = list(
    background = "#111111", panel = "#1B1B1B", grid = "#343434",
    text_main = "#F5F5F5", text_secondary = "#BFBFBF",
    accent = "#E50914",
    sequential = c("#1B1B1B", "#E50914"),
    diverging  = c("#3B82F6", "#1B1B1B", "#E50914"),
    ordinal    = c("#2ECC71", "#82E0AA", "#F1C40F", "#E67E22", "#E50914")
  ),
  github_dark = list(
    background = "#0D1117", panel = "#161B22", grid = "#30363D",
    text_main = "#C9D1D9", text_secondary = "#8B949E",
    accent = "#F85149",
    sequential = c("#161B22", "#F85149"),
    diverging  = c("#58A6FF", "#161B22", "#F85149"),
    ordinal    = c("#3FB950", "#7EE787", "#D29922", "#DB6D28", "#F85149")
  ),
  scientific = list(
    background = "#FFFFFF", panel = "#FFFFFF", grid = "#E6E6E6",
    text_main = "#161616", text_secondary = "#525252",
    accent = "#0F62FE",
    sequential = c("#FFFFFF", "#0F62FE"),
    diverging  = c("#DA1E28", "#FFFFFF", "#0F62FE"),
    ordinal    = c("#24A148", "#8DD9A8", "#F1C21B", "#EA8125", "#DA1E28")
  )
)

## fixed ordering for the categorical recovery_status factor, used everywhere
.bsst_status_levels <- c("Recovered", "Recovered (wide)", "Borderline",
                          "Not Recovered", "Severely Not Recovered")

.bsst_get_palette <- function(palette) {
  if (!palette %in% names(.bsst_palettes)) {
    stop("Unknown palette '", palette, "'. Use one of: ",
         paste(names(.bsst_palettes), collapse = ", "))
  }
  .bsst_palettes[[palette]]
}

## ============================================================
## ggplot2 theme helper
## ============================================================
.bsst_theme_ggplot <- function(palette) {
  pal <- .bsst_get_palette(palette)

  thm <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = pal$background, color = NA),
      panel.background  = ggplot2::element_rect(fill = pal$panel, color = NA),
      panel.grid.major  = ggplot2::element_line(color = pal$grid, linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_line(color = pal$grid, linewidth = 0.15),
      text              = ggplot2::element_text(color = pal$text_main),
      axis.text         = ggplot2::element_text(color = pal$text_secondary),
      axis.title        = ggplot2::element_text(color = pal$text_main),
      plot.title        = ggplot2::element_text(color = pal$text_main, face = "bold"),
      plot.subtitle     = ggplot2::element_text(color = pal$text_secondary),
      legend.background = ggplot2::element_rect(fill = pal$background, color = NA),
      legend.key        = ggplot2::element_rect(fill = pal$panel, color = NA),
      legend.text       = ggplot2::element_text(color = pal$text_secondary),
      legend.title      = ggplot2::element_text(color = pal$text_main),
      strip.background  = ggplot2::element_rect(fill = pal$panel, color = NA),
      strip.text        = ggplot2::element_text(color = pal$text_main, face = "bold")
    )

  list(theme = thm, palette = pal)
}

## ============================================================
## plotly layout helper (used by later 2D/3D/surrogate functions)
## ============================================================
.bsst_theme_plotly_layout <- function(palette) {
  pal <- .bsst_get_palette(palette)
  list(
    paper_bgcolor = pal$background,
    plot_bgcolor  = pal$panel,
    font = list(color = pal$text_main),
    xaxis = list(gridcolor = pal$grid, zerolinecolor = pal$grid, color = pal$text_secondary),
    yaxis = list(gridcolor = pal$grid, zerolinecolor = pal$grid, color = pal$text_secondary)
  )
}
## null-coalesce helper, in case not already defined elsewhere in the package
`%||%` <- function(a, b) if (is.null(a)) b else a

## ============================================================
## Internal: identify varying design variables from either
## bsst_stress() or bsst_bo() output
## ============================================================
.bsst_get_design_vars <- function(x) {
  if (!is.null(x$design_used)) {
    # bsst_stress() output
    dv <- x$design_used
    candidate_cols <- setdiff(names(dv), "design_point_id")
    varying <- candidate_cols[vapply(candidate_cols, function(nm) length(unique(dv[[nm]])) > 1, logical(1))]
    all_vars <- candidate_cols
  } else if (!is.null(x$bounds)) {
    # bsst_bo() output
    all_vars <- names(x$bounds)
    varying <- all_vars  # all declared bounds are, by construction, varying in BO
  } else {
    stop("x must be the output of bsst_stress() or bsst_bo().")
  }
  list(varying = varying, all = all_vars)
}
