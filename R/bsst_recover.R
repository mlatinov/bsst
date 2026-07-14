#' Single Parameter Recovery Test
#'
#' Runs one simulation-and-fit cycle for a Stan model and reports whether
#' each parameter of interest was recovered, based on where the true
#' (data-generating) value falls relative to the posterior distribution.
#'
#' @details
#' \strong{What this function does and does not measure.}
#' With a single simulated dataset, the quantity \code{estimate - true_value}
#' is a single realization of estimation \emph{error}, not an estimate of
#' \emph{bias}. Bias is a property of the estimator, defined as
#' \deqn{\text{Bias}(\hat\theta) = E[\hat\theta] - \theta_0,}
#' which can only be approximated by averaging errors across many
#' independent simulated datasets generated from the same \eqn{\theta_0}
#' (see \code{\link{bsst_stress}} for the repeated-simulation extension).
#' This function is a fast, single-draw diagnostic intended to catch gross
#' implementation errors (sign flips, indexing mistakes, non-identifiability,
#' severe miscalibration) before investing in a full simulation study — not
#' to certify the estimator as unbiased.
#'
#' \strong{Recovery classification.}
#' For each parameter, equal-tailed posterior credible intervals are computed
#' at each level in \code{ci_levels} (default 50\%, 80\%, 95\%, 99\%). The
#' narrowest interval that still contains the true value determines the
#' classification:
#' \itemize{
#'   \item true value inside the 50\% interval — \emph{Recovered}
#'   \item inside 80\% but not 50\% — \emph{Recovered (wide)}
#'   \item inside 95\% but not 80\% — \emph{Borderline}
#'   \item inside 99\% but not 95\% — \emph{Not Recovered}
#'   \item outside the 99\% interval — \emph{Severely Not Recovered}
#' }
#' This is nonparametric (it does not assume posterior normality) and
#' distinguishes a true value that is just outside the 95\% interval from
#' one that is off by an order of magnitude, which a single hard 95\%
#' in/out cutoff would not.
#'
#' A standardized error is also reported for continuous severity ranking:
#' \deqn{z = \frac{\bar\theta_{\text{post}} - \theta_0}{\text{sd}(\theta_{\text{post}})}}
#' where \eqn{\bar\theta_{\text{post}}} is the posterior mean and
#' \eqn{\text{sd}(\theta_{\text{post}})} the posterior standard deviation
#' of the draws for that parameter.
#'
#' \strong{Convergence gating.} Recovery results are only meaningful if the
#' sampler actually converged. Diagnostics (max R-hat, min bulk/tail ESS,
#' total divergences) are computed and compared against
#' \code{diagnostics_thresholds}; a warning is issued (not an error) if any
#' threshold is violated, since the table can still be informative for
#' debugging even when convergence failed.
#'
#' @param sim_fn An R function that simulates one dataset. Must return a
#'   list with elements \code{data} (a list or data frame passed directly to
#'   \code{cmdstanr}'s \code{$sample(data = ...)}) and \code{true_values}
#'   (a named list of the parameter values used to generate the data — names
#'   should match the corresponding Stan \code{parameters} block variable
#'   names wherever possible, to enable auto-matching).
#' @param sim_args Named list of arguments passed to \code{sim_fn}.
#' @param stan_model Either a file path to a \code{.stan} file (compiled
#'   internally via \code{cmdstanr::cmdstan_model()}) or an already-compiled
#'   \code{CmdStanModel} object (recommended when calling repeatedly, to
#'   avoid recompiling).
#' @param pars_of_interest Character vector of Stan parameter names whose
#'   recovery should be assessed and reported. Declaring this explicitly
#'   (rather than reporting all fitted parameters) avoids the "Table 2
#'   fallacy" — reporting and implicitly ascribing causal/inferential
#'   meaning to every coefficient in a model regardless of whether it was
#'   the actual target of interest.
#' @param pars_manual_map Optional named list overriding or supplementing
#'   automatic name matching between \code{true_values} and the fitted
#'   Stan variable names, e.g. \code{list(alpha = 1.5, beta = -0.8)}. Useful
#'   when a single named true value doesn't map cleanly onto the model's
#'   parameter structure (e.g., a vector \code{true_values$beta} in the sim
#'   function that actually corresponds to two separate scalar Stan
#'   parameters).
#' @param seed Optional single integer seed, passed to \code{sim_fn} (if it
#'   accepts a \code{seed} argument) and to \code{CmdStanModel$sample()}.
#' @param ci_levels Numeric vector of credible interval widths used for
#'   nested-interval recovery classification. Default
#'   \code{c(0.50, 0.80, 0.95, 0.99)}.
#' @param diagnostics_thresholds Named list of fixed convergence thresholds:
#'   \code{rhat_max}, \code{ess_bulk_min}, \code{ess_tail_min},
#'   \code{divergences_max}.
#' @param ... Additional arguments passed through to
#'   \code{CmdStanModel$sample()} (e.g. \code{chains}, \code{iter_warmup},
#'   \code{iter_sampling}, \code{parallel_chains}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{summary_table}{Data frame, one row per parameter, with columns
#'     \code{parameter}, \code{true_value}, \code{estimate_mean},
#'     \code{ci_lower}, \code{ci_upper}, \code{ci_level_used},
#'     \code{std_error}, \code{error}, \code{recovery_status}.}
#'   \item{diagnostics}{List of convergence diagnostics and a logical
#'     \code{converged} flag.}
#'   \item{fit}{The raw \code{CmdStanMCMC} fit object, for follow-up plots.}
#'   \item{data_used}{The simulated data actually passed to Stan.}
#'   \item{true_values}{The true values used, as returned by \code{sim_fn}.}
#' }
#'
#' @examples
#' \dontrun{
#' sim_linear_reg <- function(n, beta, sigma, seed = NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   x <- rnorm(n)
#'   y <- beta[1] + beta[2] * x + rnorm(n, sd = sigma)
#'   list(
#'     data = list(N = n, x = x, y = y),
#'     true_values = list(beta = beta, sigma = sigma)
#'   )
#' }
#'
#' result <- bsst_recover(
#'   sim_fn = sim_linear_reg,
#'   sim_args = list(n = 200, beta = c(1.5, -0.8), sigma = 2),
#'   stan_model = "linear_reg.stan",
#'   pars_of_interest = c("alpha", "beta", "sigma"),
#'   pars_manual_map = list(alpha = 1.5, beta = -0.8),
#'   seed = 123,
#'   chains = 4, iter_warmup = 1000, iter_sampling = 1000
#' )
#'
#' result$summary_table
#' result$diagnostics$converged
#' }
#'
#' @export
bsst_recover <- function(sim_fn,
                          sim_args = list(),
                          stan_model,
                          pars_of_interest,
                          pars_manual_map = NULL,
                          seed = NULL,
                          ci_levels = c(0.50, 0.80, 0.95, 0.99),
                          diagnostics_thresholds = list(
                            rhat_max = 1.01,
                            ess_bulk_min = 400,
                            ess_tail_min = 400,
                            divergences_max = 0
                          ),
                          ...) {

  ## ---- 1. Run simulation ----
  sim_formals <- names(formals(sim_fn))
  if (!is.null(seed)) {
    set.seed(seed)
    if ("seed" %in% sim_formals) sim_args$seed <- seed
  }
  sim_out <- do.call(sim_fn, sim_args)

  if (!all(c("data", "true_values") %in% names(sim_out))) {
    stop("sim_fn must return a list with elements 'data' and 'true_values'.")
  }

  data_used   <- sim_out$data
  true_values <- sim_out$true_values

  ## ---- 2. Compile or reuse model ----
  if (is.character(stan_model)) {
    model <- cmdstanr::cmdstan_model(stan_model)
  } else {
    model <- stan_model  # assume already a CmdStanModel
  }

  ## ---- 3. Sample ----
  fit <- model$sample(data = data_used, seed = seed, ...)

  ## ---- 4. Diagnostics ----
  fit_summary <- fit$summary()
  fit_summary <- fit_summary[fit_summary$variable != "lp__", ]

  diag_raw <- fit$diagnostic_summary(quiet = TRUE)
  total_divergences <- sum(diag_raw$num_divergent)

  diagnostics <- list(
    rhat_max     = max(fit_summary$rhat, na.rm = TRUE),
    ess_bulk_min = min(fit_summary$ess_bulk, na.rm = TRUE),
    ess_tail_min = min(fit_summary$ess_tail, na.rm = TRUE),
    divergences  = total_divergences
  )
  diagnostics$converged <- (
    diagnostics$rhat_max     <= diagnostics_thresholds$rhat_max &&
    diagnostics$ess_bulk_min >= diagnostics_thresholds$ess_bulk_min &&
    diagnostics$ess_tail_min >= diagnostics_thresholds$ess_tail_min &&
    diagnostics$divergences  <= diagnostics_thresholds$divergences_max
  )

  if (!diagnostics$converged) {
    warning("Convergence thresholds not met — recovery table may be unreliable. Check `diagnostics`.")
  }

  ## ---- 5. Match true_values to Stan parameter names ----
  draws <- fit$draws(format = "draws_df")
  draws_vars <- posterior::variables(fit$draws())

  strip_index <- function(x) sub("\\[.*\\]$", "", x)

  flat_true <- list()  # named list: stan_par_name -> true scalar value
  unmatched <- character(0)

  for (nm in names(true_values)) {
    tv <- true_values[[nm]]
    k <- length(tv)

    if (k == 1) {
      if (nm %in% draws_vars) {
        flat_true[[nm]] <- tv
      } else {
        unmatched <- c(unmatched, nm)
      }
    } else {
      expected <- paste0(nm, "[", seq_len(k), "]")
      if (all(expected %in% draws_vars)) {
        for (i in seq_len(k)) flat_true[[ expected[i] ]] <- tv[i]
      } else {
        unmatched <- c(unmatched, nm)
      }
    }
  }

  if (length(unmatched) > 0) {
    warning(paste0(
      "Could not auto-match: ", paste(unmatched, collapse = ", "),
      ". Supply these via pars_manual_map if needed."
    ))
  }

  # manual map supplements/overrides auto-matching
  if (!is.null(pars_manual_map)) {
    for (nm in names(pars_manual_map)) {
      flat_true[[nm]] <- pars_manual_map[[nm]]
    }
  }

  # keep only entries whose base name is in pars_of_interest
  keep <- vapply(names(flat_true), function(nm) strip_index(nm) %in% pars_of_interest, logical(1))
  flat_true <- flat_true[keep]

  if (length(flat_true) == 0) {
    stop("No parameters in pars_of_interest could be matched to the fitted model's variables.")
  }

  ## ---- 6. Classification helper ----
  classify_recovery <- function(true_val, draws_vec, ci_levels) {
    ci_levels <- sort(ci_levels)
    mean_est <- mean(draws_vec)
    sd_est   <- sd(draws_vec)
    std_error <- (mean_est - true_val) / sd_est
    error <- mean_est - true_val

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
      status <- switch(
        as.character(result$ci_level_used),
        "0.5"  = "Recovered",
        "0.8"  = "Recovered (wide)",
        "0.95" = "Borderline",
        "0.99" = "Not Recovered",
        paste0("Contained within ", result$ci_level_used * 100, "% CI")
      )
    }

    list(
      estimate_mean   = mean_est,
      ci_lower        = result$ci_lower,
      ci_upper        = result$ci_upper,
      ci_level_used   = result$ci_level_used,
      std_error       = std_error,
      error           = error,
      recovery_status = status
    )
  }

  ## ---- 7. Build summary table ----
  rows <- lapply(names(flat_true), function(nm) {
    true_val <- flat_true[[nm]]
    draws_vec <- draws[[nm]]
    res <- classify_recovery(true_val, draws_vec, ci_levels)

    data.frame(
      parameter       = nm,
      true_value      = true_val,
      estimate_mean   = res$estimate_mean,
      ci_lower        = res$ci_lower,
      ci_upper        = res$ci_upper,
      ci_level_used   = res$ci_level_used,
      std_error       = res$std_error,
      error           = res$error,
      recovery_status = res$recovery_status,
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, rows)
  rownames(summary_table) <- NULL

  ## ---- 8. Return ----
  list(
    summary_table = summary_table,
    diagnostics   = diagnostics,
    fit           = fit,
    data_used     = data_used,
    true_values   = true_values
  )
}