#' Simulation-Based Stress Testing Across a Design Space
#'
#' Repeatedly simulates data and fits a Stan model across a user-defined
#' grid, Latin Hypercube design, or custom set of simulation settings (e.g.
#' varying sample size, effect size, measurement error), producing a tidy
#' table describing where model recovery degrades across that design space.
#'
#' @details
#' \strong{One replicate per design point.} For computational tractability
#' (avoiding the need for cloud-scale infrastructure), this function fits
#' exactly one simulated dataset per design point rather than averaging over
#' repeated replicates at each point. Consequently the \code{error} and
#' \code{std_error} columns are single-draw estimation errors, not bias
#' estimates in the formal sense (see \code{\link{bsst_recover}} for the
#' same distinction) — a large error at a given point may reflect either a
#' genuine model weakness or ordinary sampling variability in that one
#' dataset. This function characterizes an \emph{error surface}; only
#' repeated-replicate designs (not currently implemented) would
#' characterize a true bias surface.
#'
#' \strong{Design space construction.} \code{design} is a named list where
#' each element is either a fixed scalar (held constant across all points)
#' or a vector of values to explore. Three construction modes are
#' supported via \code{design_type}:
#' \itemize{
#'   \item \code{"full_factorial"}: the full Cartesian product of all
#'     varying elements (via \code{expand.grid}). Grows multiplicatively
#'     with the number of variables and their resolution — impractical
#'     beyond 2–3 varying dimensions at fine resolution.
#'   \item \code{"lhs"}: a Latin Hypercube sample of size
#'     \code{n_design_points} over the continuous range (min, max) of each
#'     varying variable, giving good space-filling coverage at a
#'     user-controlled, fixed computational budget regardless of
#'     dimensionality.
#'   \item \code{"custom"}: a user-supplied data frame of design points
#'     (\code{design_grid}), used as-is with no modification — for
#'     externally generated designs (e.g., a specific sample-size study
#'     plan).
#' }
#'
#' \strong{Per-point recovery metrics.} For each design point and each
#' declared parameter in \code{pars_of_interest}, the same nested-interval
#' recovery classification and standardized error described in
#' \code{\link{bsst_recover}} are computed. Sampler diagnostics (R-hat,
#' bulk/tail ESS, divergences) are also recorded per point, and points that
#' fail convergence thresholds are flagged (\code{converged = FALSE}) but
#' still retained in \code{raw_table}, since a systematic pattern of
#' convergence failure across the design space is itself an important
#' finding.
#'
#' \strong{Stress score.} An optional scalar per-point objective is computed
#' from the standardized errors of only the declared \code{pars_of_interest}
#' — deliberately excluding any other model parameters, to avoid the
#' "Table 2 fallacy" of implicitly treating every fitted coefficient as
#' equally worth optimizing over. Built-in presets:
#' \deqn{\text{max\_abs\_zscore} = \max_j |z_j|}
#' \deqn{\text{mean\_sq\_zscore} = \frac{1}{J}\sum_{j=1}^J z_j^2}
#' where \eqn{z_j} is the standardized error of the \eqn{j}-th tracked
#' parameter at that design point, and \eqn{J} is the number of tracked
#' parameters. \code{max_abs_zscore} implements worst-case logic (a point
#' is flagged as bad if \emph{any} tracked parameter is poorly recovered,
#' even if others are fine); \code{mean_sq_zscore} reflects aggregate
#' severity across all tracked parameters. A user-supplied function taking
#' a numeric vector of standardized errors and returning a scalar may be
#' passed instead.
#'
#' \strong{Failure handling.} Simulation or sampling errors at a given
#' design point (e.g., a sampler crash at an extreme parameter combination)
#' are caught, logged in \code{failures_table} with the offending design
#' point and error message, and do not halt the overall run.
#'
#' @param sim_fn An R simulation function; see \code{\link{bsst_recover}}
#'   for the required return contract (\code{list(data, true_values)}).
#' @param sim_args_fixed Named list of \code{sim_fn} arguments held constant
#'   across every design point (nuisance/background parameters not being
#'   varied in this study).
#' @param design Named list specifying the design space: each element is
#'   either a fixed scalar or a vector/range of values to explore for that
#'   simulation variable (e.g. \code{list(n = c(50, 100, 200), effect_size =
#'   c(0.2, 0.5, 0.8))}).
#' @param design_type One of \code{"full_factorial"}, \code{"lhs"},
#'   \code{"custom"}.
#' @param design_grid Data frame of design points, required when
#'   \code{design_type = "custom"}; used verbatim.
#' @param n_design_points Number of points to draw, required when
#'   \code{design_type = "lhs"}.
#' @param stan_model File path or precompiled \code{CmdStanModel}.
#' @param pars_of_interest Character vector of Stan parameter names tracked
#'   for recovery and stress-score computation (see Table 2 fallacy note
#'   above).
#' @param pars_manual_map Optional manual true-value mapping; see
#'   \code{\link{bsst_recover}}.
#' @param seed Optional integer seed; each design point is assigned a
#'   deterministic derived seed (\code{seed + design_point_id}) for
#'   reproducibility across sequential and parallel execution.
#' @param ci_levels Nested credible interval widths for recovery
#'   classification.
#' @param diagnostics_thresholds Fixed convergence thresholds; see
#'   \code{\link{bsst_recover}}.
#' @param keep_fits One of \code{"none"} or \code{"all"}; whether to retain
#'   fitted \code{CmdStanMCMC} objects for every design point. Given the
#'   number of fits involved, \code{"none"} (default) is strongly
#'   recommended unless the design space is small.
#' @param parallel Logical; whether to distribute design points across
#'   parallel workers via \code{future.apply}.
#' @param n_workers Number of parallel workers, if \code{parallel = TRUE}.
#' @param objective_fn Either \code{"max_abs_zscore"}, \code{"mean_sq_zscore"},
#'   a user-supplied function, or \code{NULL} to skip stress-score
#'   computation entirely.
#' @param ... Passed through to \code{CmdStanModel$sample()}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{raw_table}{Data frame, one row per (design point, parameter),
#'     including simulation settings, recovery metrics, diagnostics, and
#'     \code{stress_score} (broadcast across all parameter rows of the same
#'     design point).}
#'   \item{failures_table}{Data frame logging design points that errored,
#'     with the offending settings and error message.}
#'   \item{design_used}{Data frame of the actual design points explored.}
#'   \item{fits}{Named list of retained fit objects, or \code{NULL} if
#'     \code{keep_fits = "none"}.}
#' }
#'
#' @seealso \code{\link{bsst_recover}} for the single-dataset diagnostic
#'   this function builds on; \code{\link{bsst_bo}} for adaptive
#'   (Bayesian-optimization-driven) search over the same design space.
#'
#' @examples
#' \dontrun{
#' sim_linear_reg <- function(n, beta, sigma, seed = NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   x <- rnorm(n)
#'   y <- beta[1] + beta[2] * x + rnorm(n, sd = sigma)
#'   list(data = list(N = n, x = x, y = y),
#'        true_values = list(beta = beta, sigma = sigma))
#' }
#'
#' stress_out <- bsst_stress(
#'   sim_fn = sim_linear_reg,
#'   sim_args_fixed = list(beta = c(1.5, -0.8)),
#'   design = list(n = c(20, 50, 100, 300), sigma = c(0.5, 1, 2, 5)),
#'   design_type = "full_factorial",
#'   stan_model = "linear_reg.stan",
#'   pars_of_interest = c("alpha", "beta", "sigma"),
#'   pars_manual_map = list(alpha = 1.5, beta = -0.8),
#'   objective_fn = "max_abs_zscore",
#'   seed = 42,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#'
#' head(stress_out$raw_table)
#' stress_out$failures_table
#' }
#'
#' @export
bsst_stress <- function(sim_fn,
                         sim_args_fixed = list(),
                         design,
                         design_type = c("full_factorial", "lhs", "custom"),
                         design_grid = NULL,
                         n_design_points = NULL,
                         stan_model,
                         pars_of_interest,
                         pars_manual_map = NULL,
                         seed = NULL,
                         ci_levels = c(0.50, 0.80, 0.95, 0.99),
                         diagnostics_thresholds = list(
                           rhat_max = 1.01, ess_bulk_min = 400,
                           ess_tail_min = 400, divergences_max = 0
                         ),
                         keep_fits = c("none", "failures", "all"),
                         parallel = FALSE,
                         n_workers = NULL,
                         objective_fn = "max_abs_zscore",
                         ...) {

  design_type <- match.arg(design_type)
  keep_fits   <- match.arg(keep_fits)

  ## resolve objective function
  obj_fn <- if (is.character(objective_fn)) {
    if (!objective_fn %in% names(.bsst_objective_presets)) {
      stop("Unknown objective_fn preset. Use one of: ",
           paste(names(.bsst_objective_presets), collapse = ", "), ", or supply a function.")
    }
    .bsst_objective_presets[[objective_fn]]
  } else if (is.function(objective_fn)) {
    objective_fn
  } else if (is.null(objective_fn)) {
    NULL
  } else stop("objective_fn must be a string, function, or NULL.")

  ## compile model once
  model <- if (is.character(stan_model)) cmdstanr::cmdstan_model(stan_model) else stan_model

  ## build design points
  design_used <- .bsst_generate_design(design, design_type, design_grid, n_design_points)
  varying_names <- setdiff(names(design_used), "design_point_id")
  # only the truly-varying columns (drop constants introduced from `fixed`) for point extraction
  # (fixed ones just get folded into sim_args_fixed at call time, harmless either way)

  point_ids <- design_used$design_point_id

  run_one <- function(pid) {
    row <- design_used[design_used$design_point_id == pid, , drop = FALSE]
    point_values <- as.list(row[ , setdiff(names(row), "design_point_id"), drop = FALSE])

    pt_seed <- if (!is.null(seed)) seed + pid else NULL  # deterministic, non-colliding per point

    res <- .bsst_simulate_and_fit_point(
      design_point_values = point_values,
      sim_fn = sim_fn, sim_args_fixed = sim_args_fixed,
      model = model, pars_of_interest = pars_of_interest,
      pars_manual_map = pars_manual_map, seed = pt_seed,
      ci_levels = ci_levels, diagnostics_thresholds = diagnostics_thresholds,
      keep_fit = (keep_fits == "all"), ...
    )

    if (!res$success && keep_fits == "failures") {
      # nothing to keep, fit never completed; just ensure failure gets logged below
    }

    list(design_point_id = pid, result = res, point_values = point_values)
  }

  ## ---- execute (sequential or parallel) ----
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' required for parallel = TRUE.")
    }
    if (!is.null(n_workers)) future::plan(future::multisession, workers = n_workers)
    all_results <- future.apply::future_lapply(point_ids, run_one, future.seed = TRUE)
  } else {
    all_results <- lapply(point_ids, run_one)
  }

  ## ---- assemble outputs ----
  raw_rows <- list()
  failure_rows <- list()
  kept_fits <- list()

  for (r in all_results) {
    pid <- r$design_point_id
    res <- r$result

    if (!res$success) {
      failure_rows[[length(failure_rows) + 1]] <- data.frame(
        design_point_id = pid,
        as.data.frame(r$point_values, stringsAsFactors = FALSE),
        error_msg = res$error_msg,
        stringsAsFactors = FALSE
      )
      next
    }

    par_rows <- res$par_rows
    par_rows$design_point_id <- pid

    ## objective: computed ONLY over pars_of_interest std_error, per point
    if (!is.null(obj_fn)) {
      par_rows$stress_score <- obj_fn(par_rows$std_error)
    }

    raw_rows[[length(raw_rows) + 1]] <- cbind(
      as.data.frame(r$point_values, stringsAsFactors = FALSE),
      par_rows
    )

    if (keep_fits == "all" && !is.null(res$fit)) kept_fits[[as.character(pid)]] <- res$fit
  }

  raw_table <- if (length(raw_rows) > 0) do.call(rbind, raw_rows) else data.frame()
  failures_table <- if (length(failure_rows) > 0) do.call(rbind, failure_rows) else data.frame()

  list(
    raw_table      = raw_table,
    failures_table = failures_table,
    design_used    = design_used,
    fits           = if (length(kept_fits) > 0) kept_fits else NULL
  )
}