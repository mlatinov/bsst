// Data Input Layer 
data{
    // Observations Prior switch and Indexes
    int <lower = 0, upper = 1> prior_only;
    int <lower = 2 > N_lakes;
    int <lower = 2>  N_weeks;
    int <lower = 1>  dt;

    // Covariates 
    vector[N_lakes] flow_rate;
    vector[N_lakes] initial_spill_size;

    // Outcome
    matrix[N_weeks, N_lakes] pollution; 
}

// Data transformations 
transformed data {
   // Center the initial spill size 
   vector[N_lakes] c_initial_spill_size;
   c_initial_spill_size = initial_spill_size - mean(initial_spill_size); 

}
// Model Parameters 
parameters{
    // Initial State Parameters 
    real mean_polution_concentration;
    real beta_spill_size;
    real <lower = 0.001 > sd_initial_state;
    vector[N_lakes] z_intitial;

    // Rate of Change Paramters 
    real mean_rate_decay;
    real beta_flow_rate;
    real <lower = 0.001> sd_rate_decay;
    vector[N_lakes] z_rate;

    // Observational Paramters for the OU Process 
    real <lower = 0.001> sd_obs;
    real <lower = 0.001> sd_ou;
    matrix[N_weeks, N_lakes] C_raw;

}
// Transform Model Parameters
transformed parameters {
    // Combine and create the SDE paramters for initial state and rate of change 
    vector[N_lakes] po;
    vector[N_lakes] r;
    for(i in 1:N_lakes){
        po[i] = mean_polution_concentration + beta_spill_size * c_initial_spill_size[i] + sd_initial_state * z_intitial[i];
        r[i]  = exp(mean_rate_decay + beta_flow_rate * flow_rate[i] + sd_rate_decay * z_rate[i]);
    }
    // Latent OU Process
    matrix[N_weeks, N_lakes] C;
    // Initialize the process 
    C[1, ] = to_row_vector(po);
    // OU Process 
    for(i in 2:N_weeks){
        vector[N_lakes] mean_t = to_vector(C[i - 1, ]) .* exp(-r * dt);
        vector[N_lakes] sd_t   = sqrt((sd_ou^2 ./ (2 * r)) .* (1 - exp(-2 * r * dt)));
        C[i ,] = to_row_vector(mean_t + sd_t .* to_vector(C_raw[i, ]));
    }
}
// Model 
model{
    // Initial State Priors 
    mean_polution_concentration ~ normal(10, 2);
    beta_spill_size ~ normal(0.8, 0.5);
    sd_initial_state ~ exponential(1);
    z_intitial ~ std_normal();

    // Rate of Change Priors
    mean_rate_decay ~ normal(0.08, 0.01);
    beta_flow_rate  ~ normal(0.03, 0.01);
    sd_rate_decay ~ exponential(1);
    z_rate ~ std_normal();

    // Observational Priors for the OU Process 
    sd_obs ~ exponential(1);
    sd_ou  ~ exponential(1);
    to_vector(C_raw) ~ std_normal();

    // Observation Model
    if(prior_only == 0){
        to_vector(pollution) ~ normal(to_vector(C), sd_obs);
    } 
}
// Minimal Generated Quantites 
generated quantities{
    matrix[N_weeks, N_lakes] pollution_rep;
    matrix[N_weeks, N_lakes] log_lik;
    for(i in 1:N_weeks){
        for(j in 1:N_lakes){
            pollution_rep[i, j] = normal_rng(C[i, j], sd_obs);
            log_lik[i, j] = normal_lpdf(pollution[i, j] | C[i, j], sd_obs);
        }
    }
}