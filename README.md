# BSST — Bayesian Simulation-based Severe Testing

**BSST** is an R package for validating Stan models through simulation-based
parameter recovery experiments. It answers a question every applied
Bayesian modeler eventually has to face, usually the hard way:

> *If I already knew the true parameter values, would my model actually
> recover them?*

BSST doesn't just check whether a model *can* recover its parameters on a
single convenient dataset. It's built to help you find the conditions
under which it **can't** — the sample sizes, effect sizes, group counts,
and design choices where inference quietly breaks down — before those
conditions show up in real data.

---

## Table of Contents

- [Why BSST](#why-bsst)
- [Core Concepts](#core-concepts)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [The Simulation Function Contract](#the-simulation-function-contract)
- [Core Functions](#core-functions)
  - [`bsst_recover()`](#bsst_recover)
  - [`bsst_stress()`](#bsst_stress)
  - [`bsst_bo()`](#bsst_bo)
- [Plotting Functions](#plotting-functions)
- [Palettes](#palettes)
- [Package Philosophy & Design Constraints](#package-philosophy--design-constraints)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)

---

## Why BSST

Fitting a Stan model and getting clean diagnostics (R-hat, ESS, no
divergences) tells you the sampler did its job. It tells you **nothing**
about whether the model itself can recover the truth — that's a separate,
harder question, and one that's rarely tested systematically.

BSST is built around three escalating layers of that question:

1. **Does it work at all, once?** — a fast, single-simulation sanity check.
2. **Where does it break, across a design space?** — a systematic sweep
   over simulation conditions (sample size, effect size, group counts,
   measurement error, anything your simulation exposes).
3. **Where is the single worst case?** — an adaptive search that spends
   its compute budget hunting for the region of the design space where
   recovery is most severely wrong, rather than exhaustively gridding
   everything.

Each layer builds on the one before it, and all three share the same
underlying simulate-fit-classify machinery, the same recovery
classification logic, and the same plotting API.

---

## Core Concepts

A few ideas run through the entire package. Understanding them up front
will save you from misreading the output.

### Error vs. Bias

With a single simulated dataset, `estimate − true_value` is one
realization of estimation **error**, not an estimate of **bias**. Bias is
a property of an estimator — `E[θ̂] − θ₀` — and can only be approximated by
averaging errors across many independent replicates at the same true
value.

BSST is explicit about this distinction throughout:

- `bsst_recover()` and `bsst_stress()` report **error** and **standardized
  error**, computed from a single replicate per point, for cost reasons.
- Nothing in the package currently averages across replicates, so nothing
  in the package currently estimates bias in the formal sense. A single
  large error at one point in a stress sweep may reflect a genuine model
  weakness, or it may just be sampling noise in that one dataset — the
  tables and plots characterize an **error surface**, not a bias surface.

If you need true bias estimates, you'll need to average multiple
replicates per condition yourself (or wait for a future BSST function
built around that).

### The Table 2 Fallacy

Reporting recovery statistics for *every* fitted model parameter,
regardless of which ones actually matter for your research question, is a
well-known trap ("Table 2 fallacy") that implicitly treats every
coefficient as equally deserving of scrutiny and interpretation.

BSST forces you to declare `pars_of_interest` explicitly in every core
function. Objective functions used for stress scoring and Bayesian
optimization are computed **only** over these declared parameters — a
poorly-recovered nuisance parameter can never silently drive the
optimizer's search or inflate a stress score you weren't asking about.

### Nested Credible Interval Classification

Rather than a single hard 95%-in/95%-out cutoff, BSST classifies recovery
using nested credible intervals (default: 50/80/95/99%). The narrowest
interval containing the true value determines the label:

| True value falls within... | Label |
|---|---|
| 50% interval | Recovered |
| 80% but not 50% | Recovered (wide) |
| 95% but not 80% | Borderline |
| 99% but not 95% | Not Recovered |
| Outside 99% | Severely Not Recovered |

This is nonparametric (no normality assumption) and distinguishes a true
value that's just outside the 95% interval from one that's off by an
order of magnitude — information a single cutoff would throw away.

### Convergence Gates Everything

Recovery results from a fit that didn't converge are not meaningful. Every
core function computes standard sampler diagnostics (max R-hat, min
bulk/tail ESS, divergences) against fixed thresholds and flags
non-convergence — but does not discard the fit or hide the numbers,
since a systematic pattern of convergence failure across a design space
is itself an important finding (often the *most* important one).

**Always check `diagnostics$converged` — or the `converged` column in a
stress/BO table — before interpreting any recovery number.**

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("mlatinov/bsst")
```

BSST requires a working CmdStan installation via `cmdstanr`:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

### Dependencies

| Package | Role | Required? |
|---|---|---|
| `cmdstanr` | Stan interface | Hard dependency |
| `posterior` | Draws handling | Hard dependency |
| `ggplot2` | Static plotting | Hard dependency |
| `DiceKriging` | GP surrogate for `bsst_bo()` | Hard dependency |
| `lhs` | Latin Hypercube designs | Hard dependency |
| `plotly` | Interactive/3D plotting | Suggested — only needed for `engine = "plotly"` and 3D functions |
| `future.apply` | Parallel execution | Suggested — only needed for `parallel = TRUE` |

---

## Quick Start

```r
library(bsst)

# 1. Write a simulation function returning list(data, true_values)
sim_linear_reg <- function(n, beta, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- rnorm(n)
  y <- beta[1] + beta[2] * x + rnorm(n, sd = sigma)
  list(
    data = list(N = n, x = x, y = y),
    true_values = list(beta = beta, sigma = sigma)
  )
}

# 2. Compile your Stan model once
model <- cmdstanr::cmdstan_model("linear_reg.stan")

# 3. Check recovery on a single simulated dataset
result <- bsst_recover(
  sim_fn = sim_linear_reg,
  sim_args = list(n = 200, beta = c(1.5, -0.8), sigma = 2),
  stan_model = model,
  pars_of_interest = c("alpha", "beta", "sigma"),
  pars_manual_map = list(alpha = 1.5, beta = -0.8),
  seed = 123,
  chains = 4, iter_warmup = 1000, iter_sampling = 1000
)

bsst_plot_recovery(result)
bsst_plot_diagnostics(result)
```

---

## The Simulation Function Contract

Every BSST function is built around one convention: **you supply an R
function that simulates one dataset and returns the ground truth used to
generate it.**

```r
sim_fn <- function(...) {
  # ... simulate a dataset ...
  list(
    data        = list(...),   # passed directly to CmdStanModel$sample(data = ...)
    true_values = list(...)    # named list; names should match Stan `parameters` block names
  )
}
```

- **`data`** must already be shaped exactly as your Stan model's `data`
  block expects (index arrays, counts, everything). BSST does no
  reshaping of its own — if your underlying simulator (e.g. an external
  package) returns a flat data frame, write a thin adapter that reshapes
  it, derives any counts (`N`, group counts) from the actual simulated
  data rather than hardcoding them, and renames columns to match your
  Stan variable names.
- **`true_values`** should use names that match your Stan `parameters`
  block wherever possible, so BSST can auto-match them to posterior draws.
  Vector-valued parameters are given as a single named element containing
  a vector (e.g. `true_values$beta <- c(0.5, -1.2)`), and are auto-expanded
  against `beta[1], beta[2], ...` in the fitted draws.
- When a true value can't be auto-matched (mismatched structure, a vector
  in simulation that maps to several scalar Stan parameters, etc.), supply
  `pars_manual_map` to fill or override the mapping explicitly. BSST warns
  rather than guesses when a match fails.
- If your simulator's underlying function takes many arguments you don't
  want to hardcode into the adapter, just forward `...`:
  `function(...) { raw <- external_sim(...); list(data = ..., true_values = ...) }`.
  This also makes the same adapter reusable across `bsst_recover()`,
  `bsst_stress()`, and `bsst_bo()` without modification.

---

## Core Functions

### `bsst_recover()`

**Question answered:** *Did this one fit work?*

Runs a single simulate-and-fit cycle and reports, per declared parameter,
the posterior mean, the narrowest containing credible interval, the
standardized error, and a qualitative recovery label — gated by
convergence diagnostics.

```r
result <- bsst_recover(
  sim_fn, sim_args = list(),
  stan_model,                 # file path or precompiled CmdStanModel
  pars_of_interest,
  pars_manual_map = NULL,
  seed = NULL,
  ci_levels = c(0.50, 0.80, 0.95, 0.99),
  diagnostics_thresholds = list(rhat_max = 1.01, ess_bulk_min = 400,
                                 ess_tail_min = 400, divergences_max = 0),
  ...                          # forwarded to CmdStanModel$sample()
)

result$summary_table   # one row per tracked parameter
result$diagnostics     # rhat/ess/divergences + converged flag
result$fit             # raw CmdStanMCMC, for further inspection
```

### `bsst_stress()`

**Question answered:** *Across a range of simulation conditions, where
does recovery degrade?*

Repeats the single-fit cycle across a grid, a Latin Hypercube sample, or a
user-supplied custom design of simulation settings — varying sample size,
effect size, group counts, measurement error, or anything else your
simulation function exposes as an argument.

```r
stress_out <- bsst_stress(
  sim_fn, sim_args_fixed = list(),
  design = list(n = c(50, 100, 200), effect_size = c(0.2, 0.5, 0.8)),
  design_type = c("full_factorial", "lhs", "custom"),
  stan_model, pars_of_interest, pars_manual_map = NULL,
  objective_fn = "max_abs_zscore",   # or "mean_sq_zscore", or a custom function
  keep_fits = "none",
  parallel = FALSE, n_workers = NULL,
  ...
)

stress_out$raw_table       # one row per (design point, parameter)
stress_out$failures_table  # logged errors, run doesn't halt on a bad point
stress_out$design_used     # the actual points explored
```

`objective_fn` collapses the standardized errors of your declared
`pars_of_interest` into a single per-point `stress_score`:

- `"max_abs_zscore"` (default) — worst-case across tracked parameters.
- `"mean_sq_zscore"` — aggregate severity across tracked parameters.
- Or supply your own `function(std_error_vec) -> scalar`.

### `bsst_bo()`

**Question answered:** *Where, in a continuous design space, is recovery
worst — without exhaustively gridding it?*

Uses Bayesian optimization (a noise-aware Gaussian process surrogate via
`DiceKriging`, batch proposals via the Kriging Believer heuristic, and
Expected Improvement as the acquisition function) to actively search a
continuous design space for the region of highest `stress_score`,
optionally warm-started from a `bsst_stress()` run.

```r
bo_out <- bsst_bo(
  sim_fn, sim_args_fixed = list(),
  design = list(n = c(20, 300)),     # continuous search box, min/max per variable
  integer_vars = "n",                 # rounded after each continuous proposal
  stan_model, pars_of_interest,
  warm_start = stress_out,            # optional; recommended
  n_iter = 8, batch_size = 4,
  objective_fn = "max_abs_zscore",    # must be smooth; avoid discrete objectives here
  ...
)

bo_out$best_point   # the single worst design point found
bo_out$surrogate     # fitted GP, for surface plotting
bo_out$raw_table      # every evaluated point (warm start + init + BO batches)
```

**Note:** because each stress score is computed from a single simulated
dataset per point, the surrogate is fit with an explicit nugget (noise)
term — a noise-free GP here would overfit sampling noise and chase
spurious spikes rather than genuine high-error regions.

---

## Plotting Functions

All plotting functions are prefixed `bsst_plot_*`, share a common
`palette` argument, and — where a 2D/3D view is offered — a common
`engine = "ggplot2" | "plotly"` argument.

| Function | Consumes | Answers |
|---|---|---|
| `bsst_plot_recovery()` | `bsst_recover()` | Which parameters missed, and by how much? |
| `bsst_plot_diagnostics()` | `bsst_recover()` | Did the sampler actually converge? |
| `bsst_plot_stress_1d()` | `bsst_stress()` / `bsst_bo()` | As I vary one setting, where does recovery break down? |
| `bsst_plot_stress_2d()` | `bsst_stress()` / `bsst_bo()` | Is there an interaction between two conditions? |
| `bsst_plot_stress_3d()` | `bsst_stress()` / `bsst_bo()` | Same as above, with height as a third channel (plotly only) |
| `bsst_plot_profile()` | `bsst_stress()` / `bsst_bo()` | 1D view, with an option to average rather than fix other dimensions |
| `bsst_plot_frontier()` | `bsst_stress()` / `bsst_bo()` | What's the minimum sample size / maximum tolerable error for acceptable inference? |
| `bsst_plot_status()` | `bsst_stress()` / `bsst_bo()` | Per-parameter categorical recovery status across the design space |
| `bsst_plot_diagnostic_surface()` | `bsst_stress()` / `bsst_bo()` | Is a "bad" region actually a sampler-convergence problem instead? |
| `bsst_plot_surrogate()` | `bsst_bo()` | What does the full continuous error surface look like, including unexplored regions? |
| `bsst_plot_bo_trace()` | `bsst_bo()` | Has the search converged, or does it need more iterations? |
| `bsst_plot_bo_overlay()` | `bsst_bo()` | Did the optimizer concentrate its search where it should have? |

Functions requiring exactly two design variables (`_2d`, `_3d`,
`_status`, `_surrogate`) require `fix_others` whenever the underlying
design has more than the plotted variables varying — BSST will error
with the names of what's missing rather than silently average over an
unaddressed dimension. `parameter` must always be supplied explicitly;
there is no silent "plot everything" default.

---

## Palettes

Every plot function accepts `palette = "dark_research" | "github_dark" |
"scientific"`:

- **`dark_research`** — near-black background, Netflix-red accent. Good
  for live demos and dashboards.
- **`github_dark`** — GitHub's dark theme colors. Good for docs/READMEs
  rendered on GitHub.
- **`scientific`** — white background, IBM Carbon Design System colors.
  Built for figures going into papers or reports; check your target
  venue's specific figure requirements (e.g. grayscale-safety) before
  assuming any colored palette clears them.

---

## Package Philosophy & Design Constraints

A few decisions are baked into the package on purpose and are worth
knowing about rather than fighting against:

- **One replicate per design/BO point, by default.** True bias and
  empirical coverage require averaging over many replicates per
  condition — BSST currently trades that rigor for tractable compute cost
  at scale. Treat `bsst_stress()`/`bsst_bo()` output as an *error surface*,
  not a *bias surface*, and re-run a point with `bsst_recover()` (or
  average manually) before treating any single result as conclusive.
- **Fit objects are not retained by default** (`keep_fits = "none"`) in
  `bsst_stress()`/`bsst_bo()`, since retaining hundreds of `CmdStanMCMC`
  objects is not memory-practical. Only summary statistics and
  diagnostics survive per point unless you explicitly opt in.
- **Objective functions only ever see `pars_of_interest`.** This is not
  configurable per-call in a way that lets nuisance parameters leak in —
  it's a structural guard against the Table 2 fallacy, not a default you
  can accidentally disable.
- **Integer design variables are handled by rounding, not true
  mixed-integer optimization**, in both the LHS design generator and
  `bsst_bo()`'s continuous search. This is a deliberate simplicity
  tradeoff; it can occasionally produce duplicate evaluations at coarse
  resolution.
- **Failures don't halt a run.** A simulation or sampler crash at one
  design point is caught, logged with its settings and error message, and
  the run continues — because a pattern of *where* fits fail is often as
  informative as the recovery numbers themselves.

---

## Roadmap

Known gaps, in rough priority order:

- [ ] A repeated-replicate mode (true bias / empirical coverage estimation
      across multiple simulated datasets per design point), for users
      willing to pay the added compute cost.
- [ ] Per-point origin tracking in `bsst_bo()`'s output (warm start vs.
      initial design vs. specific BO batch), currently not recorded
      internally, which limits `bsst_plot_bo_overlay()`.
- [ ] Mixed-integer-aware Bayesian optimization, replacing the current
      round-after-propose heuristic.
- [ ] A fast, small toy-model test suite, separate from slow real-model
      integration/smoke tests, for asserting exact expected behavior on
      every release.
- [ ] Formal handling of failed-point feedback into the BO surrogate
      (currently dropped entirely rather than treated as a penalized
      observation).

---

## Contributing

Issues and PRs welcome. If you're adding a new plotting function or
objective function preset, please keep to the existing conventions:
explicit `parameter`/`fix_others` arguments with no silent defaults,
shared palette/theme helpers rather than one-off styling, and clear
errors over best-effort guesses whenever a request is ambiguous.

## License

MIT — see `LICENSE`.
