test_that("BSST Core Functions Pipeline on 3L Nested Hierarchical Model", {

  # Simulate data from Rsims Library 
  adapter_simulate_schools <- function(...) {
    raw <- rsims::simulate_schools(...)
    list(
      data = list(
        N = nrow(raw),
        prior_only = 0,
        test_scores = raw$test_scores,
        ses_index   = raw$ses_index,
        class_size  = raw$class_size,
        study_h     = raw$study_hours,
        school_id   = raw$schools,
        district_id = raw$districts,
        state_id    = raw$states,
        n_school_j   = length(unique(raw$schools)),
        n_district_k = length(unique(raw$districts)),
        n_state_l    = length(unique(raw$states))
    ),
    true_values = list(
      alpha           = 70,
      sigma_state     = 2,
      sigma_district  = 3,
      sigma_school    = 5,
      ses_index_bar   = 6,
      sigma_ses_index = 1,
      beta_study_h    = 4,
      beta_class_size = -0.3,
      test_score_sigma = 6
    )
  )
}
  # Compile a Stan Model for the data 
  stan_model <- cmdstanr::cmdstan_model(stan_file = "Stan/m_hierarchical_student_achievement.stan")

  # Test the BSST Recover Function 
  result <- bsst_recover(
    sim_fn = adapter_simulate_schools,
    sim_args = list(),
    stan_model = stan_model,
    pars_of_interest = c("alpha", "sigma_state", "sigma_district", "sigma_school",
                          "ses_index_bar", "sigma_ses_index", "beta_study_h", "beta_class_size"),
    seed = 42,
    chains = 4, 
    iter_sampling = 200
  )
  # Run all the structure and Logic tests 
  expect_bsst_result(result)

  # Test BSST Stress
  # Without knowing now many districtrics we need to get stable estimands we run stress test over a 
  # grid variing the num of states LHC example
  stress_out <- bsst_stress(
    sim_fn = adapter_simulate_schools,
    sim_args_fixed = list(
      num_districts_per_state = 6,
      num_schools_per_district = 5,
      num_of_student_per_schools = 5,
      class_size_range = c(18, 32),
      mean_study_hours = 2.5,
      sd_study_hours = 0.8,
      s_beta_ses = 0.45,
      baseline_test_score = 70,
      sd_test_scores_districts = 3,
      sd_test_scores_schools = 5,
      beta_study_h = 4,
      beta_class_size = -0.3,
      beta_bar_ses_index = 6,
      sd_beta_ses_index = 1,
      test_score_sd = 6
    ),
    design = list(num_states = c(4, 8, 10, 12)),
    design_type = "full_factorial",
    stan_model = stan_model,
    pars_of_interest = c("alpha", "sigma_state", "sigma_district", "sigma_school",
                        "ses_index_bar", "sigma_ses_index", "beta_study_h", "beta_class_size"),
    seed = 1,
    objective_fn = "max_abs_zscore",
    chains = 2,
    iter_sampling = 200
  )
  # Run the tests for the structure of stress out 
  expect_recovery_grid_result(stress_out)
})