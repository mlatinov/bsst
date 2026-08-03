# tests/testthat/helper-bsst.R
expect_bsst_result <- function(result){
    expect_type(result, "list")
    expect_named(
        result,
        c(
            "summary_table",
            "diagnostics",
            "fit",
            "data_used",
            "true_values"
        )
    )
    expect_s3_class(result$summary_table, "data.frame")
    expect_s3_class(result$fit, "CmdStanMCMC")
    expect_true(
        all(result$summary_table$ci_lower <=
              result$summary_table$estimate_mean)
    )
    expect_true(
        all(result$summary_table$estimate_mean <=
              result$summary_table$ci_upper)
    )
    expect_false(any(is.na(result$summary_table)))
}
# test bsst_stress helper 
expect_recovery_grid_result <- function(result){

  # Main object structure
  testthat::expect_type(result, "list")

  testthat::expect_named(
    result,
    c(
      "raw_table",
      "failures_table",
      "design_used",
      "fits"
    )
  )


  # Raw recovery table
  testthat::expect_s3_class(
    result$raw_table,
    "data.frame"
  )

  required_columns <- c(
    "num_states",
    "parameter",
    "true_value",
    "estimate_mean",
    "error",
    "std_error",
    "ci_lower",
    "ci_upper",
    "ci_level_used",
    "recovery_status",
    "rhat",
    "ess_bulk",
    "ess_tail",
    "divergences",
    "converged",
    "design_point_id",
    "stress_score"
  )

  testthat::expect_true(
    all(required_columns %in% names(result$raw_table))
  )


  # Basic numerical sanity
  testthat::expect_true(
    all(
      result$raw_table$ci_lower <=
        result$raw_table$estimate_mean
    )
  )

  testthat::expect_true(
    all(
      result$raw_table$estimate_mean <=
        result$raw_table$ci_upper
    )
  )


  testthat::expect_false(
    anyNA(result$raw_table)
  )


  # Design grid
  testthat::expect_s3_class(
    result$design_used,
    "data.frame"
  )

  testthat::expect_true(
    all(
      c("num_states", "design_point_id") %in%
        names(result$design_used)
    )
  )


  # Failures table can legitimately be empty
  testthat::expect_s3_class(
    result$failures_table,
    "data.frame"
  )


  # Fits can be NULL because saving fits is optional
  testthat::expect_true(
    is.null(result$fits) ||
      is.list(result$fits)
  )

}
# tests/testthat/helper-bo.R
expect_bo_result <- function(result){

  # Main object structure
  testthat::expect_type(result, "list")

  testthat::expect_named(
    result,
    c(
      "raw_table",
      "failures_table",
      "best_point",
      "surrogate",
      "bounds",
      "fits"
    )
  )

  # ------------------------------------------------------------------
  # Raw results table
  # ------------------------------------------------------------------

  testthat::expect_s3_class(
    result$raw_table,
    "data.frame"
  )

  required_columns <- c(
    "num_libraries",
    "parameter",
    "true_value",
    "estimate_mean",
    "error",
    "std_error",
    "ci_lower",
    "ci_upper",
    "ci_level_used",
    "recovery_status",
    "rhat",
    "ess_bulk",
    "ess_tail",
    "divergences",
    "converged",
    "design_point_id",
    "stress_score"
  )

  testthat::expect_true(
    all(required_columns %in% names(result$raw_table))
  )

  testthat::expect_true(
    all(
      result$raw_table$ci_lower <=
        result$raw_table$estimate_mean
    )
  )

  testthat::expect_true(
    all(
      result$raw_table$estimate_mean <=
        result$raw_table$ci_upper
    )
  )

  testthat::expect_false(
    anyNA(result$raw_table)
  )

  # ------------------------------------------------------------------
  # Failures table
  # ------------------------------------------------------------------

  testthat::expect_s3_class(
    result$failures_table,
    "data.frame"
  )

  # ------------------------------------------------------------------
  # Best point
  # ------------------------------------------------------------------

  testthat::expect_s3_class(
    result$best_point,
    "data.frame"
  )

  testthat::expect_equal(
    nrow(result$best_point),
    1
  )

  testthat::expect_true(
    all(required_columns %in% names(result$best_point))
  )

  # ------------------------------------------------------------------
  # Surrogate model
  # ------------------------------------------------------------------

  testthat::expect_s4_class(
    result$surrogate,
    "km"
  )

  # ------------------------------------------------------------------
  # Bounds
  # ------------------------------------------------------------------

  testthat::expect_type(
    result$bounds,
    "list"
  )

  testthat::expect_true(
    length(result$bounds) > 0
  )

  testthat::expect_true(
    all(vapply(result$bounds, length, integer(1)) == 2)
  )

  # ------------------------------------------------------------------
  # Fits
  # ------------------------------------------------------------------

  testthat::expect_true(
    is.null(result$fits) ||
      is.list(result$fits)
  )
}