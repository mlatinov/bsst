#' Adaptive Bayesian Optimization Search for Worst-Case Model Recovery
#'
#' Actively searches a continuous simulation design space for the region(s)
#' where a Stan model's parameter recovery is worst, using Bayesian
#' optimization with a Gaussian process surrogate fit over an
#' error/stress-score surface rather than exhaustively gridding the space.
#'
#' @details
#' \strong{Why Bayesian optimization here.} \code{\link{bsst_stress}}
#' characterizes recovery over a pre-specified grid or space-filling design,
#' but its resolution is fixed in advance and computational cost grows with
#' the number of points regardless of whether most of them turn out to be
#' uninteresting. This function instead treats "where is recovery worst"
#' as a black-box maximization problem over the \code{stress_score}
#' surface and uses Bayesian optimization (BO) to concentrate expensive
#' Stan fits in regions the surrogate model believes are most promising
#' (highest expected stress), refining its belief after every batch of
#' fits.
#'
#' \strong{Surrogate model and noise.} Each evaluation of the stress score
#' at a design point comes from a \emph{single} simulated dataset (for the
#' same compute-cost reasons as \code{\link{bsst_stress}}), so repeated
#' evaluation at the same point would not return the same value — the
#' observed stress score is a noisy sample of an underlying, unobserved
#' error surface. Standard noise-free GP interpolation (as commonly used in
#' textbook BO) would overfit this sampling noise and chase spurious
#' spikes rather than genuine high-error regions. To address this, the
#' surrogate is a Gaussian process with an explicit nugget (noise) term,
#' fit via \code{DiceKriging::km(..., nugget.estim = TRUE)}, which
#' separates estimated observation noise from the smooth trend in the
#' response surface.
#'
#' \strong{Acquisition function.} Points are proposed by maximizing
#' Expected Improvement (EI) over the current surrogate, since the
#' objective here is \emph{maximization} of the stress score (larger
#' standardized error = worse recovery = more interesting):
#' \deqn{
#'   EI(x) = (\mu(x) - y^{*})\,\Phi(z) + \sigma(x)\,\phi(z), \quad
#'   z = \frac{\mu(x) - y^{*}}{\sigma(x)}
#' }
#' where \eqn{\mu(x)} and \eqn{\sigma(x)} are the surrogate's posterior
#' mean and standard deviation at candidate point \eqn{x}, \eqn{y^{*}} is
#' the best (highest) stress score observed so far, and \eqn{\Phi}, \eqn{\phi}
#' are the standard normal CDF and PDF. EI is optimized over the unit
#' hypercube via multi-start L-BFGS-B, since it can be multimodal.
#'
#' \strong{Batch proposals (Kriging Believer).} Because Stan fits are the
#' computational bottleneck and can be run in parallel, this function
#' proposes \code{batch_size} points per iteration rather than one point at
#' a time. Within a batch, after each point is chosen, its outcome is
#' temporarily "fantasized" as the surrogate's own predicted mean at that
#' point (the Kriging Believer heuristic), the GP is refit including this
#' fantasy observation, and the next point in the batch is proposed against
#' the updated surrogate — this prevents all points in a batch from
#' collapsing onto the same single maximum. Once a full batch of candidate
#' points is selected, all of them are evaluated for real (actual
#' simulation + Stan fit), and the fantasies are discarded in favor of the
#' true observed stress scores before the next iteration's surrogate fit.
#'
#' \strong{Search bounds and discreteness.} The search occurs in the
#' continuous box defined by \code{range()} of each variable declared in
#' \code{design} (matching the same declaration style as
#' \code{\link{bsst_stress}}, but interpreted as a continuous
#' \code{c(min, max)} bound rather than a discrete set of values). Variables
#' named in \code{integer_vars} (e.g. sample size) are rounded to the
#' nearest integer after being proposed in continuous space — a simple
#' approach that can occasionally cause duplicate evaluations after
#' rounding, in exchange for not requiring a mixed-integer-aware BO
#' implementation.
#'
#' \strong{Warm starting.} If \code{warm_start} is supplied (the output of
#' \code{\link{bsst_stress}}, run with a non-\code{NULL} \code{objective_fn}
#' matching the one used here), its explored design points and their
#' stress scores are used as the surrogate's initial training data,
#' avoiding redundant re-simulation of an initial space-filling design. If
#' omitted, \code{n_init} random points are drawn uniformly from the search
#' box and evaluated first.
#'
#' @param sim_fn,sim_args_fixed,stan_model,pars_of_interest,pars_manual_map,seed,ci_levels,diagnostics_thresholds
#'   As in \code{\link{bsst_stress}}.
#' @param design Named list; each element is a \code{c(min, max)} range (or
#'   any vector, from which \code{range()} is taken) defining the
#'   continuous search box for that simulation variable.
#' @param objective_fn \code{"max_abs_zscore"}, \code{"mean_sq_zscore"}, or
#'   a user function. Restricted to objectives expected to vary smoothly
#'   over the design space, since the GP surrogate requires a reasonably
#'   continuous response; discrete/step objectives (e.g. a raw coverage
#'   violation indicator) are not suitable here and should instead be
#'   examined descriptively via \code{\link{bsst_stress}}.
#' @param integer_vars Character vector naming which \code{design} variables
#'   must be rounded to the nearest integer after proposal (e.g.
#'   \code{"n"}).
#' @param warm_start Optional output of \code{\link{bsst_stress}}, used to
#'   initialize the surrogate's training data.
#' @param n_init Number of random initial points evaluated if
#'   \code{warm_start} is not supplied.
#' @param n_iter Number of BO iterations (batches), not individual points.
#' @param batch_size Number of points proposed and evaluated per iteration.
#' @param n_restarts_ei Number of multi-start restarts used when optimizing
#'   the Expected Improvement acquisition function per proposed point.
#' @param keep_fits \code{"none"} or \code{"all"}.
#' @param parallel,n_workers As in \code{\link{bsst_stress}}; particularly
#'   relevant here since each batch's real evaluations can be distributed
#'   across workers.
#' @param ... Passed through to \code{CmdStanModel$sample()}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{raw_table}{Data frame of every evaluated point (from warm start,
#'     initial design, and all BO batches), same row structure as
#'     \code{\link{bsst_stress}}'s \code{raw_table} — stackable/comparable
#'     with it.}
#'   \item{failures_table}{Logged simulation/fitting failures.}
#'   \item{best_point}{The single design point with the highest observed
#'     \code{stress_score} found during the search.}
#'   \item{surrogate}{The final fitted \code{DiceKriging} GP model, for
#'     follow-up surface plotting.}
#'   \item{bounds}{The search box actually used, per variable.}
#'   \item{fits}{Retained fit objects, if \code{keep_fits = "all"}.}
#' }
#'
#' @seealso \code{\link{bsst_stress}} for the non-adaptive grid/LHS
#'   counterpart this function can warm-start from.
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
#' # optional: warm start from a coarse grid first
#' grid_out <- bsst_stress(
#'   sim_fn = sim_linear_reg,
#'   sim_args_fixed = list(beta = c(1.5, -0.8)),
#'   design = list(n = c(20, 300), sigma = c(0.5, 5)),
#'   design_type = "lhs", n_design_points = 15,
#'   stan_model = "linear_reg.stan",
#'   pars_of_interest = c("alpha", "beta", "sigma"),
#'   pars_manual_map = list(alpha = 1.5, beta = -0.8),
#'   seed = 1
#' )
#'
#' bo_out <- bsst_bo(
#'   sim_fn = sim_linear_reg,
#'   sim_args_fixed = list(beta = c(1.5, -0.8)),
#'   design = list(n = c(20, 300), sigma = c(0.5, 5)),
#'   stan_model = "linear_reg.stan",
#'   pars_of_interest = c("alpha", "beta", "sigma"),
#'   pars_manual_map = list(alpha = 1.5, beta = -0.8),
#'   integer_vars = "n",
#'   warm_start = grid_out,
#'   n_iter = 8, batch_size = 4,
#'   seed = 2, chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#'
#' bo_out$best_point
#' }
#'
#' @export
bsst_bo <- function(sim_fn,
                     sim_args_fixed = list(),
                     design,                        # named list, element = c(min, max) per varying var
                     stan_model,
                     pars_of_interest,
                     pars_manual_map = NULL,
                     seed = NULL,
                     ci_levels = c(0.50, 0.80, 0.95, 0.99),
                     diagnostics_thresholds = list(
                       rhat_max = 1.01, ess_bulk_min = 400,
                       ess_tail_min = 400, divergences_max = 0
                     ),
                     objective_fn = "max_abs_zscore",   # restricted to continuous presets or user fn
                     integer_vars = character(0),
                     warm_start = NULL,               # optional: output of bsst_stress()
                     n_init = 10,                      # used only if warm_start is NULL
                     n_iter = 10,                       # number of BATCHES, not individual points
                     batch_size = 4,
                     n_restarts_ei = 20,
                     keep_fits = c("none", "all"),
                     parallel = FALSE,
                     n_workers = NULL,
                     ...) {

  keep_fits <- match.arg(keep_fits)

  ## ---- resolve objective function (must be continuous-friendly; coverage_violation excluded) ----
  obj_fn <- if (is.character(objective_fn)) {
    presets <- list(
      max_abs_zscore = function(v) max(abs(v), na.rm = TRUE),
      mean_sq_zscore = function(v) mean(v^2, na.rm = TRUE)
    )
    if (!objective_fn %in% names(presets)) {
      stop("objective_fn preset must be 'max_abs_zscore' or 'mean_sq_zscore' for BO (needs a smooth surface), or supply your own function.")
    }
    presets[[objective_fn]]
  } else if (is.function(objective_fn)) objective_fn
  else stop("objective_fn must be a string preset or a function.")

  ## ---- bounds from design ----
  bounds <- lapply(design, function(v) c(min(v), max(v)))
  varying_names <- names(bounds)
  d <- length(varying_names)
  if (d == 0) stop("design must declare at least one variable with a c(min, max) range.")

  fixed_extra <- list()  # nothing extra here; sim_args_fixed handles true fixed nuisance params

  ## ---- compile model once ----
  model <- if (is.character(stan_model)) cmdstanr::cmdstan_model(stan_model) else stan_model

  ## ---- helper: evaluate one real point (calls the same internal fitter as bsst_stress) ----
  next_id_counter <- new.env()
  next_id_counter$id <- 0

  eval_point <- function(point_values) {
    next_id_counter$id <- next_id_counter$id + 1
    pid <- next_id_counter$id
    pt_seed <- if (!is.null(seed)) seed + pid else NULL

    res <- .bsst_simulate_and_fit_point(
      design_point_values = point_values,
      sim_fn = sim_fn, sim_args_fixed = sim_args_fixed,
      model = model, pars_of_interest = pars_of_interest,
      pars_manual_map = pars_manual_map, seed = pt_seed,
      ci_levels = ci_levels, diagnostics_thresholds = diagnostics_thresholds,
      keep_fit = (keep_fits == "all"), ...
    )

    if (!res$success) {
      return(list(design_point_id = pid, point_values = point_values,
                  success = FALSE, error_msg = res$error_msg,
                  stress_score = NA_real_, par_rows = NULL, fit = NULL))
    }

    par_rows <- res$par_rows
    par_rows$design_point_id <- pid
    score <- obj_fn(par_rows$std_error)
    par_rows$stress_score <- score

    list(design_point_id = pid, point_values = point_values, success = TRUE,
         error_msg = NA_character_, stress_score = score,
         par_rows = cbind(as.data.frame(point_values, stringsAsFactors = FALSE), par_rows),
         fit = res$fit)
  }

  ## ---- initialize training data: warm start or fresh random/LHS design ----
  history_rows <- list()
  failure_rows <- list()
  kept_fits <- list()
  X_unit <- matrix(nrow = 0, ncol = d)
  y_obs  <- numeric(0)

  add_to_history <- function(res) {
    if (!res$success) {
      failure_rows[[length(failure_rows) + 1]] <<- data.frame(
        design_point_id = res$design_point_id,
        as.data.frame(res$point_values, stringsAsFactors = FALSE),
        error_msg = res$error_msg, stringsAsFactors = FALSE)
      return(invisible(NULL))
    }
    history_rows[[length(history_rows) + 1]] <<- res$par_rows
    if (keep_fits == "all" && !is.null(res$fit)) {
      kept_fits[[as.character(res$design_point_id)]] <<- res$fit
    }
    x_row <- .bsst_normalize(res$point_values, bounds)
    X_unit <<- rbind(X_unit, x_row)
    y_obs  <<- c(y_obs, res$stress_score)
  }

  if (!is.null(warm_start)) {
    ## pull unique design points + their stress_score from bsst_stress() output
    rt <- warm_start$raw_table
    if (is.null(rt) || !"stress_score" %in% names(rt)) {
      stop("warm_start must be the output of bsst_stress() with a non-NULL objective_fn.")
    }
    unique_pts <- rt[!duplicated(rt$design_point_id), c(varying_names, "stress_score"), drop = FALSE]
    for (i in seq_len(nrow(unique_pts))) {
      pv <- as.list(unique_pts[i, varying_names, drop = FALSE])
      x_row <- .bsst_normalize(pv, bounds)
      X_unit <- rbind(X_unit, x_row)
      y_obs  <- c(y_obs, unique_pts$stress_score[i])
    }
  } else {
    ## fresh random initial design (simple uniform sampling within bounds)
    init_points <- lapply(seq_len(n_init), function(i) {
      pv <- .bsst_denormalize(runif(d), bounds, integer_vars)
      pv
    })
    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      if (!is.null(n_workers)) future::plan(future::multisession, workers = n_workers)
      init_results <- future.apply::future_lapply(init_points, eval_point, future.seed = TRUE)
    } else {
      init_results <- lapply(init_points, eval_point)
    }
    for (r in init_results) add_to_history(r)
  }

  if (nrow(X_unit) < 2) stop("Need at least 2 successfully evaluated initial points to fit a surrogate. Increase n_init or check for simulation/fit failures.")

  ## ============================================================
  ## Main BO loop — batch proposals via Kriging Believer
  ## ============================================================
  for (iter in seq_len(n_iter)) {

    km_model <- tryCatch(
      DiceKriging::km(design = as.data.frame(X_unit), response = y_obs,
                       nugget.estim = TRUE, control = list(trace = FALSE)),
      error = function(e) {
        warning("km() fit failed at iteration ", iter, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(km_model)) break

    y_best <- max(y_obs)

    ## ---- batch proposal: Kriging Believer ----
    batch_x_unit <- matrix(nrow = 0, ncol = d)
    fantasy_X <- X_unit
    fantasy_y <- y_obs
    fantasy_model <- km_model

    for (b in seq_len(batch_size)) {
      x_next <- .bsst_propose_next(fantasy_model, max(fantasy_y), d, n_restarts = n_restarts_ei)
      batch_x_unit <- rbind(batch_x_unit, x_next)

      if (b < batch_size) {
        pred <- predict(fantasy_model, newdata = as.data.frame(matrix(x_next, nrow = 1)),
                         type = "UK", checkNames = FALSE)
        fantasy_X <- rbind(fantasy_X, x_next)
        fantasy_y <- c(fantasy_y, pred$mean)  # believer: assume predicted mean is the outcome
        fantasy_model <- tryCatch(
          DiceKriging::km(design = as.data.frame(fantasy_X), response = fantasy_y,
                           nugget.estim = TRUE, control = list(trace = FALSE)),
          error = function(e) fantasy_model  # if refit fails, reuse previous surrogate
        )
      }
    }

    ## ---- evaluate the real batch (actual Stan fits) ----
    batch_points <- lapply(seq_len(nrow(batch_x_unit)), function(i) {
      .bsst_denormalize(batch_x_unit[i, ], bounds, integer_vars)
    })

    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      if (!is.null(n_workers)) future::plan(future::multisession, workers = n_workers)
      batch_results <- future.apply::future_lapply(batch_points, eval_point, future.seed = TRUE)
    } else {
      batch_results <- lapply(batch_points, eval_point)
    }

    for (r in batch_results) add_to_history(r)
  }

  ## ============================================================
  ## Assemble outputs
  ## ============================================================
  raw_table <- if (length(history_rows) > 0) do.call(rbind, history_rows) else data.frame()
  failures_table <- if (length(failure_rows) > 0) do.call(rbind, failure_rows) else data.frame()

  best_row <- if (nrow(raw_table) > 0) {
    unique_pts <- raw_table[!duplicated(raw_table$design_point_id), ]
    unique_pts[which.max(unique_pts$stress_score), ]
  } else NULL

  final_surrogate <- tryCatch(
    DiceKriging::km(design = as.data.frame(X_unit), response = y_obs,
                     nugget.estim = TRUE, control = list(trace = FALSE)),
    error = function(e) NULL
  )

  list(
    raw_table       = raw_table,
    failures_table  = failures_table,
    best_point      = best_row,
    surrogate       = final_surrogate,
    bounds          = bounds,
    fits            = if (length(kept_fits) > 0) kept_fits else NULL
  )
}