## ============================================================
## Internal: normalize / denormalize design points to [0,1]^d
## ============================================================
.bsst_normalize <- function(x_row, bounds) {
  vapply(names(bounds), function(nm) {
    (x_row[[nm]] - bounds[[nm]][1]) / (bounds[[nm]][2] - bounds[[nm]][1])
  }, numeric(1))
}

.bsst_denormalize <- function(x_unit, bounds, integer_vars) {
  vals <- lapply(seq_along(bounds), function(i) {
    nm <- names(bounds)[i]
    v <- x_unit[i] * (bounds[[nm]][2] - bounds[[nm]][1]) + bounds[[nm]][1]
    if (nm %in% integer_vars) v <- round(v)
    v
  })
  names(vals) <- names(bounds)
  vals
}

## ============================================================
## Internal: Expected Improvement (for MAXIMIZATION, since higher
## stress_score = worse recovery = what we're hunting for)
## ============================================================
.bsst_ei <- function(x_unit_mat, km_model, y_best, epsilon = 0) {
  pred <- predict(km_model, newdata = as.data.frame(x_unit_mat), type = "UK", checkNames = FALSE)
  mu    <- pred$mean
  sigma <- pred$sd
  sigma[sigma < 1e-8] <- 1e-8

  z <- (mu - y_best - epsilon) / sigma
  ei <- (mu - y_best - epsilon) * pnorm(z) + sigma * dnorm(z)
  ei[sigma < 1e-7] <- 0
  ei
}

## ============================================================
## Internal: optimize EI over unit hypercube via multi-start L-BFGS-B
## ============================================================
.bsst_propose_next <- function(km_model, y_best, d, n_restarts = 20) {
  best_val <- -Inf
  best_x <- NULL

  starts <- matrix(runif(n_restarts * d), nrow = n_restarts, ncol = d)

  for (i in seq_len(n_restarts)) {
    obj <- function(x) -.bsst_ei(matrix(x, nrow = 1), km_model, y_best)
    res <- tryCatch(
      optim(starts[i, ], obj, method = "L-BFGS-B",
            lower = rep(0, d), upper = rep(1, d)),
      error = function(e) NULL
    )
    if (!is.null(res) && -res$value > best_val) {
      best_val <- -res$value
      best_x <- res$par
    }
  }

  if (is.null(best_x)) best_x <- runif(d)  # fallback: random point if all optim calls failed
  best_x
}
## ============================================================
## Internal: recovery classification (shared logic w/ bsst_recover)
## ============================================================
.bsst_classify_recovery <- function(true_val, draws_vec, ci_levels) {
  ci_levels <- sort(ci_levels)
  mean_est  <- mean(draws_vec)
  sd_est    <- sd(draws_vec)

  std_error <- if (sd_est > 0) (mean_est - true_val) / sd_est else NA_real_
  error     <- mean_est - true_val

  result <- NULL
  for (level in ci_levels) {
    alpha <- 1 - level
    lo <- unname(quantile(draws_vec, probs = alpha / 2))
    hi <- unname(quantile(draws_vec, probs = 1 - alpha / 2))
    if (true_val >= lo && true_val <= hi) {
      result <- list(ci_lower = lo, ci_upper = hi, ci_level_used = level)
      break
    }
  }

  if (is.null(result)) {
    widest <- max(ci_levels)
    alpha <- 1 - widest
    lo <- unname(quantile(draws_vec, probs = alpha / 2))
    hi <- unname(quantile(draws_vec, probs = 1 - alpha / 2))
    result <- list(ci_lower = lo, ci_upper = hi, ci_level_used = widest)
    status <- "Severely Not Recovered"
  } else {
    status <- switch(as.character(result$ci_level_used),
      "0.5" = "Recovered", "0.8" = "Recovered (wide)",
      "0.95" = "Borderline", "0.99" = "Not Recovered",
      paste0("Contained within ", result$ci_level_used * 100, "% CI"))
  }

  list(estimate_mean = mean_est, ci_lower = result$ci_lower,
       ci_upper = result$ci_upper, ci_level_used = result$ci_level_used,
       std_error = std_error, error = error, recovery_status = status)
}

## ============================================================
## Internal: match true_values to draws (same logic as bsst_recover)
## ============================================================
.bsst_match_true_values <- function(true_values, draws_vars, pars_of_interest, pars_manual_map) {
  strip_index <- function(x) sub("\\[.*\\]$", "", x)
  flat_true <- list()
  unmatched <- character(0)

  for (nm in names(true_values)) {
    tv <- true_values[[nm]]
    k <- length(tv)
    if (k == 1) {
      if (nm %in% draws_vars) flat_true[[nm]] <- tv else unmatched <- c(unmatched, nm)
    } else {
      expected <- paste0(nm, "[", seq_len(k), "]")
      if (all(expected %in% draws_vars)) {
        for (i in seq_len(k)) flat_true[[ expected[i] ]] <- tv[i]
      } else {
        unmatched <- c(unmatched, nm)
      }
    }
  }

  if (!is.null(pars_manual_map)) {
    for (nm in names(pars_manual_map)) flat_true[[nm]] <- pars_manual_map[[nm]]
  }

  keep <- vapply(names(flat_true), function(nm) strip_index(nm) %in% pars_of_interest, logical(1))
  list(flat_true = flat_true[keep], unmatched = unmatched)
}

## ============================================================
## Internal: simulate + fit ONE design point, return lean stats only
## This is the function BOTH bsst_stress() and the future BO function will call.
## ============================================================
.bsst_simulate_and_fit_point <- function(design_point_values,   # named list, the varying args for this point
                                          sim_fn, sim_args_fixed,
                                          model, pars_of_interest, pars_manual_map,
                                          seed, ci_levels, diagnostics_thresholds,
                                          keep_fit = FALSE, ...) {

  sim_args <- modifyList(sim_args_fixed, design_point_values)

  sim_formals <- names(formals(sim_fn))
  if (!is.null(seed)) {
    set.seed(seed)
    if ("seed" %in% sim_formals) sim_args$seed <- seed
  }

  out <- tryCatch({

    sim_out <- do.call(sim_fn, sim_args)
    if (!all(c("data", "true_values") %in% names(sim_out))) {
      stop("sim_fn must return list(data = ..., true_values = ...)")
    }
    data_used   <- sim_out$data
    true_values <- sim_out$true_values

    fit <- model$sample(data = data_used, seed = seed, refresh = 0, ...)

    fit_summary <- fit$summary()
    fit_summary <- fit_summary[fit_summary$variable != "lp__", ]
    diag_raw <- fit$diagnostic_summary(quiet = TRUE)

    diagnostics <- list(
      rhat_max     = max(fit_summary$rhat, na.rm = TRUE),
      ess_bulk_min = min(fit_summary$ess_bulk, na.rm = TRUE),
      ess_tail_min = min(fit_summary$ess_tail, na.rm = TRUE),
      divergences  = sum(diag_raw$num_divergent)
    )
    diagnostics$converged <- (
      diagnostics$rhat_max     <= diagnostics_thresholds$rhat_max &&
      diagnostics$ess_bulk_min >= diagnostics_thresholds$ess_bulk_min &&
      diagnostics$ess_tail_min >= diagnostics_thresholds$ess_tail_min &&
      diagnostics$divergences  <= diagnostics_thresholds$divergences_max
    )

    draws <- fit$draws(format = "draws_df")
    draws_vars <- posterior::variables(fit$draws())

    matched <- .bsst_match_true_values(true_values, draws_vars, pars_of_interest, pars_manual_map)
    flat_true <- matched$flat_true

    if (length(flat_true) == 0) stop("No pars_of_interest matched fitted variables.")

    par_rows <- lapply(names(flat_true), function(nm) {
      res <- .bsst_classify_recovery(flat_true[[nm]], draws[[nm]], ci_levels)
      data.frame(parameter = nm, true_value = flat_true[[nm]],
                 estimate_mean = res$estimate_mean, error = res$error,
                 std_error = res$std_error, ci_lower = res$ci_lower,
                 ci_upper = res$ci_upper, ci_level_used = res$ci_level_used,
                 recovery_status = res$recovery_status,
                 rhat = diagnostics$rhat_max, ess_bulk = diagnostics$ess_bulk_min,
                 ess_tail = diagnostics$ess_tail_min, divergences = diagnostics$divergences,
                 converged = diagnostics$converged,
                 stringsAsFactors = FALSE)
    })

    list(success = TRUE, par_rows = do.call(rbind, par_rows),
         fit = if (keep_fit) fit else NULL, error_msg = NA_character_)

  }, error = function(e) {
    list(success = FALSE, par_rows = NULL, fit = NULL, error_msg = conditionMessage(e))
  })

  out
}

## ============================================================
## Internal: design generation
## ============================================================
.bsst_generate_design <- function(design, design_type, design_grid, n_design_points) {

  varying <- design[vapply(design, length, integer(1)) > 1]
  fixed   <- design[vapply(design, length, integer(1)) == 1]

  if (design_type == "custom") {
    if (is.null(design_grid)) stop("design_grid required when design_type = 'custom'.")
    grid <- design_grid

  } else if (design_type == "full_factorial") {
    if (length(varying) == 0) {
      grid <- data.frame(dummy = 1)[, FALSE, drop = FALSE]
      grid[1, ] <- NA  # single point, all fixed
      grid <- as.data.frame(matrix(nrow = 1, ncol = 0))
    } else {
      grid <- expand.grid(varying, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }

  } else if (design_type == "lhs") {
    if (length(varying) == 0) stop("design_type = 'lhs' requires at least one varying design variable.")
    if (is.null(n_design_points)) stop("n_design_points required when design_type = 'lhs'.")

    raw_lhs <- lhs::randomLHS(n_design_points, length(varying))
    grid <- as.data.frame(raw_lhs)
    names(grid) <- names(varying)

    # scale unit hypercube [0,1] to each variable's declared range
    for (nm in names(varying)) {
      rng <- range(varying[[nm]])
      grid[[nm]] <- grid[[nm]] * (rng[2] - rng[1]) + rng[1]
      # if original values looked like integers (e.g., sample sizes), round
      if (all(varying[[nm]] == round(varying[[nm]]))) grid[[nm]] <- round(grid[[nm]])
    }

  } else {
    stop("design_type must be one of 'full_factorial', 'lhs', 'custom'.")
  }

  # attach fixed values as constant columns
  for (nm in names(fixed)) grid[[nm]] <- fixed[[nm]]

  grid$design_point_id <- seq_len(nrow(grid))
  grid
}

## ============================================================
## Default objective functions (operate on std_error vector for
## pars_of_interest ONLY, at a single design point)
## ============================================================
.bsst_objective_presets <- list(
  max_abs_zscore    = function(std_error_vec) max(abs(std_error_vec), na.rm = TRUE),
  mean_sq_zscore    = function(std_error_vec) mean(std_error_vec^2, na.rm = TRUE)
)
