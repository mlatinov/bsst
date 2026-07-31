// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 2> N_soil_types;
    int <lower = 0, upper = 1> prior_only;
    array[N] int soil_type_id;

    // Splines
    int <lower = 2> N_basis; 
    matrix[N, N_basis] B_fertilizer;

    // Outcome 
    vector[N] biomass;
}

// Model Parameters
parameters{
    // Random Intercept Parameters
    real baseline_soil_type_biomass;
    real <lower = 0.001> soil_type_biomass_sd;
    vector[N_soil_types] zj;

    // Fertilezer P-spline Parameters Population 
    real wp1;
    real wp2;
    real <lower = 0.001> sd_wp;
    vector[N_basis - 2] zp;

    // Fertilier P-spline Parameters per unit diviations
    matrix[N_basis, N_soil_types] z_sp;
    real <lower = 0.001> omega;

    // Observation model parameters
    real <lower = 0.001> sd_obs;

}
// Transform Parameters
transformed parameters {
   // Ensemble the random intercept 
   vector[N_soil_types] alpha_j;
   for(i in 1:N_soil_types){
    alpha_j[i] = baseline_soil_type_biomass + soil_type_biomass_sd * zj[i];
   }

   // Ensemble the P-spline 
   vector[N_basis] w_pop;
   w_pop[1] = wp1;
   w_pop[2] = wp2;
   for(k in 3:N_basis){
    w_pop[k] = 2 * w_pop[k -1] - w_pop[k - 2] + sd_wp * zp[k- 2];
   }

   // NCP
   matrix[N_basis, N_soil_types] w_unit;
   for(i in 1:N_soil_types){
    w_unit[,i] = w_pop + omega .* z_sp[, i];
   }
}
// Model 
model{
    // Random Intercept Priors
    baseline_soil_type_biomass ~ normal(0, 1);
    soil_type_biomass_sd ~ exponential(1);
    zj ~ std_normal();

    // Fertilezer P-spline Priors Population 
    wp1 ~ normal(0, 1);
    wp2 ~ normal(0, 1);
    sd_wp ~ exponential(1);
    zp    ~ std_normal();

    // Fertilier P-spline Priors per unit diviations
    to_vector(z_sp) ~ std_normal();
    omega ~ exponential(1);

    // Observation model priors
    sd_obs ~ exponential(1);

    // Model Likelihood
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = alpha_j[soil_type_id[i]] + B_fertilizer[i] * w_unit[, soil_type_id[i]];
        }
    // Observation Model 
    biomass ~ normal(mu, sd_obs);
    }
}
// Minimal Generated Quantities
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i] = alpha_j[soil_type_id[i]] + B_fertilizer[i] * w_unit[, soil_type_id[i]];
        log_lik[i] = normal_lpdf(biomass[i] | mu[i], sd_obs);
        y_rep[i]   = normal_rng(mu[i], sd_obs);
    }
}