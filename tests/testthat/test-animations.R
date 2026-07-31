
# Fit the model 
data  <- simulate_ode_bacterial_growth()
model <- cmdstanr::cmdstan_model(stan_file = "Stan/m_ode_bacterial_growth.stan") 
culture_lookup <- data[!duplicated(data$culture_id), ]
culture_lookup <- culture_lookup[order(culture_lookup$culture_id), ]

fit <- model$sample(
  data = list(
    N = nrow(data),
    prior_only = 0,
    N_cultures = length(unique(data$culture_id)),
    cultures_id = data$culture_id,
    day_id      = data$day,
    temp_c      = culture_lookup$temperature,
    nutrient_c  = culture_lookup$nutrient_c,
    bacterial_count = data$bacterial_count
  ),
  seed = 42,
  chains = 4,
  iter_sampling = 4000
)

draws_df <- fit$draws("y_rep") |> posterior::as_draws_df()
draws_thin <- posterior::subset_draws(draws_df, draw = sample(unique(draws_df$.draw), 1500))
row_lookup <- data.frame(
  row  = seq_len(nrow(data)),
  time = data$day,
  unit = data$culture_id
)

long <- draws_thin |>
  dplyr::select(-.chain, -.iteration) |>
  dplyr::rename(draw = .draw) |>
  tidyr::pivot_longer(
    starts_with("y_rep"),
    names_to = "row",
    values_to = "value"
  ) |>
  dplyr::mutate(
    row = as.integer(gsub("\\D", "", row))
  ) |>
  dplyr::left_join(row_lookup, by = "row")

animate_temporal_posterior(long, unit = 2, mode = "cloud", theme = "midnight",
                            observed = data.frame(time = data$day[data$culture_id==1],
                                                   value = data$bacterial_count[data$culture_id==1]))
