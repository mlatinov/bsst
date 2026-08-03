test_that("BSST Core Functions Pipeline on 1L Nested Hierarchical Model", {

  ## Simulate data from Rsims library 
  adapter <- function(...){
    raw <- rsims::simulate_libraries(...)
    list(
      data = list(
        N          = nrow(raw),
        prior_only = 0,
        is_school_holiday = raw$is_school_holiday,
        rain_mm           = raw$rain_mm,
        held_book         = raw$held_book,
        staff_count       = raw$staff_count,
        n_libraries       = length(unique(raw$library_id)),
        library_id        = raw$library_id, 
        weekly_book_loans = raw$weekly_book_loans
      ),
      true_values = list(
        mean_book_loans  = log(140),
        sigma_book_count = 0.2,
        c_beta_holiday   = 0.18,
        c_beta_rain      = 0.005,
        c_beta_held      = 0.04,
        c_beta_staff     = 0.06
      )
    )
  }

  ## Compile Stan model for the data 
  stan_model <- cmdstanr::cmdstan_model("Stan/hierarchical_models/m_hierarchical_libraries.stan")

  # Test the BSST Recover Function 
  result <- bsst_recover(
    sim_fn = adapter,
    sim_args = list(),
    stan_model = stan_model,
    pars_of_interest = c("mean_book_loans", "sigma_book_count", "c_beta_holiday", "c_beta_rain","c_beta_held", "c_beta_staff"),
    seed = 42,
    chains = 4, 
    iter_sampling = 2000,
    adapt_delta   = 0.95
  )
  # Run all the structure and Logic tests 
  expect_bsst_result(result)

  # Test BSST Stress without knowing knowing now many libraries w
  stress_out <- bsst_stress(
    sim_fn = adapter,
    sim_args_fixed = list(
       num_libraries = 20,
       num_weeks = 52,
       rain_mm_rate = 0.08,
       baseline_held_book = log(8),
       staff_count_ranges = c(3, 10),
       h_beta_holiday = 0.35,
       h_beta_rain = 0.015,
       mean_book_loans = log(140),
       sigma_book_count = 0.2,
       c_beta_holiday = 0.18,
       c_beta_rain = 0.005,
       c_beta_held = 0.04,
       c_beta_staff = 0.06
    ),
    design = list(num_libraries = c(2, 4, 8, 16)),
    design_type = "full_factorial",
    stan_model = stan_model,
    pars_of_interest = c("mean_book_loans", "sigma_book_count", "c_beta_holiday", "c_beta_rain","c_beta_held", "c_beta_staff"),
    seed = 42,
    objective_fn = "max_abs_zscore",
    chains = 4,
    iter_sampling = 2000
  )
  # Run the tests for the structure of stress out 
  expect_recovery_grid_result(stress_out)

  # Test BSST BO 
  bo_out <- bsst_bo(
    sim_fn = adapter,
    sim_args_fixed = list(
      num_weeks = 52,
      rain_mm_rate = 0.08,
      baseline_held_book = log(8),
      staff_count_ranges = c(3, 10),
      h_beta_holiday = 0.35,
      h_beta_rain = 0.015,
      mean_book_loans = log(140),
      sigma_book_count = 0.2,
      c_beta_holiday = 0.18,
      c_beta_rain = 0.005,
      c_beta_held = 0.04,
      c_beta_staff = 0.06
    ),
    design = list(num_libraries = c(2, 16)),   # continuous search box bounds
    stan_model = stan_model,
    pars_of_interest = c("mean_book_loans", "sigma_book_count", "c_beta_holiday",
                          "c_beta_rain", "c_beta_held", "c_beta_staff"),
    integer_vars = "num_libraries",
    warm_start = stress_out,
    seed = 42,
    objective_fn = "max_abs_zscore",
    n_iter = 3,
    batch_size = 2,
    chains = 4,
    iter_sampling = 2000,
    adapt_delta = 0.95
  )

  ## Test the output of BO 
  expect_bo_result(bo_out)  

})