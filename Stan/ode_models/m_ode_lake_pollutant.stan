// Data Input Layer
data{
    // Observation and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Indexing 
    int <lower = 2> N_lakes;
    array[N] int    lake_id;
    array[N] int    week_id;
    
    // Covariates
    vector[N_lakes] flow;
    vector[N_lakes] intial_spill;
    
    // Outcome 
    vector[N] polution;
}
// Model Parameters 
parameters{
    // ODE parameters //

    // Initial State parameters 
    real u_polution_c;
    real <lower = 0.001> sd_initial_state;
    vector[N_lakes] z_polution;

    // Rate of Change parameters
    real u_rate;
    real <lower = 0.001> sd_rate;
    vector[N_lakes] z_rate;

    // Covariate parameters 
    real beta_flow;
    real beta_spill;

    // Outcome Observation Paramters 
    real <lower = 0.001> sd_obs;
}
// Transform Paramters to recover the the ODE coefs
transformed parameters {
   vector[N_lakes] log_xo; // Initial State
   vector[N_lakes] log_r;  // Rate of Change
   vector[N_lakes] xo;
   vector[N_lakes] r; 
   for(i in 1:N_lakes){
    // Initial State X0 per Lake 
    log_xo[i] = u_polution_c + beta_spill * log(intial_spill[i]) +  sd_initial_state * z_polution[i];
    // Rate of change R per Lake 
    log_r[i]  = u_rate + beta_flow * flow[i] + sd_rate * z_rate[i];
   }
   xo = exp(log_xo);
   r  = exp(log_r);
}
// Model 
model{
    // Initial State Priors 
    u_polution_c ~ normal(log(10), 0.5);
    sd_initial_state ~ exponential(1);
    z_polution       ~ std_normal();

    // Rate of Change Priors
    u_rate ~ normal(log(0.08), 0.5);
    sd_rate ~ exponential(1);
    z_rate  ~ std_normal();

    // Covariate parameters 
    beta_flow ~ normal(0.03, 0.1);
    beta_spill ~ normal(0.8, 1);

    // Outcome Observation Paramters 
    sd_obs ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            // Analytical Solution to the ODE 
            mu[i] = xo[lake_id[i]] * exp(-r[lake_id[i]] * week_id[i]);
        }
    // Sample the observations from LogNormal distribution 
    polution ~ lognormal(log(mu), sd_obs);
    }
}
// Minimal Generated Quantities
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i]      = xo[lake_id[i]] * exp(-r[lake_id[i]] * week_id[i]);
        log_lik[i] = lognormal_lpdf(polution[i] | log(mu[i]), sd_obs);
        y_rep[i]   = lognormal_rng(log(mu[i]), sd_obs);
    }
}