// Data Layer Input 
data{
    // Observations prior switch and Indexing 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;
    int <lower = 2> N_subjects;
    array[N] int hour_id;
    array[N] int subject_id;

    // Covariates 
    vector[N_subjects] ambient_temp;
    vector[N_subjects] activity_level;

    // Outcome 
    vector[N] body_temperature;
}
// Transform Data 
transformed data {
    // Center the covariates 
    vector[N_subjects] c_ambient_temp;
    vector[N_subjects] c_activity_level;
    c_ambient_temp   = ambient_temp - mean(ambient_temp);
    c_activity_level = activity_level - mean(activity_level);
}
// Model Parameters 
parameters{
    // Parameters to construct a
    real mean_heat_inflow_rate;
    real <lower = 0.002> sd_subjects_heat_inflow_rate;
    real beta_ambient;
    vector[N_subjects] z_a;

    // Parameters to construct b 
    real mean_outflow_rate;
    real <lower = 0.001> sd_subjects_outflow_rate;
    real beta_activity;
    vector[N_subjects] z_b;

    // Paramters for the initial state 
    real mean_initial_state;
    real <lower = 0.001> sd_subjects_initial_state;
    vector[N_subjects] z_xo;

    // Parameters for the Observation Model 
    real <lower = 0.001> sd_obs;
}
// Transform Model Parameters 
transformed parameters {
    // Reconstruct the ODE equation terms from the paramters 
    vector[N_subjects] log_a;
    vector[N_subjects] log_b;
    vector[N_subjects] xo;
    vector[N_subjects] a;
    vector[N_subjects] b;
    for(i in 1:N_subjects){
        log_a[i] = mean_heat_inflow_rate + beta_ambient  * c_ambient_temp[i] + sd_subjects_heat_inflow_rate * z_a[i];
        log_b[i] = mean_outflow_rate     + beta_activity * c_activity_level[i] + sd_subjects_outflow_rate * z_b[i];
        xo[i]    = mean_initial_state + sd_subjects_initial_state * z_xo[i];
    }
    a = exp(log_a);
    b = exp(log_b);
}
// Model
model{
    // Priors to construct a
    mean_heat_inflow_rate ~ normal(log(0.15) + log(37), 0.5);
    sd_subjects_heat_inflow_rate ~ exponential(1);
    beta_ambient ~ normal(0.015, 0.01);
    z_a          ~ std_normal();

    // Priors to construct b 
    mean_outflow_rate ~ normal(log(0.15), 0.5);
    sd_subjects_outflow_rate ~ exponential(1);
    beta_activity ~ normal(-0.03, 0.1);
    z_b           ~ std_normal();

    // Priors for the initial state 
    mean_initial_state ~ normal(30, 5);
    sd_subjects_initial_state ~ exponential(1);
    z_xo                      ~ std_normal();

    // Priors for the Observation Model 
    sd_obs ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = 
                (a[subject_id[i]] / b[subject_id[i]])
                + 
               (xo[subject_id[i]] - (a[subject_id[i]] / b[subject_id[i]]))
                *
                exp(-b[subject_id[i]] * hour_id[i]);
        }
    // Observational Model 
    body_temperature ~ lognormal(log(mu), sd_obs);
    }

}
// Minimal Generated Quantities 
generated quantities {
    vector[N] log_lik;
    vector[N] y_rep;
    for (i in 1:N) {
        real eq_i  = a[subject_id[i]] / b[subject_id[i]];
        real mu_i  = eq_i + (xo[subject_id[i]] - eq_i) * exp(-b[subject_id[i]] * hour_id[i]);
        log_lik[i] = lognormal_lpdf(body_temperature[i] | log(mu_i), sd_obs);
        y_rep[i]   = lognormal_rng(log(mu_i), sd_obs);
    }
}