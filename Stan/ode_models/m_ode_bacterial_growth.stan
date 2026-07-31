// Data Input Layer
data{
    // Observation N Indexing and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;
    int <lower = 2> N_cultures;
    array[N] int cultures_id;
    array[N] int day_id;

    // Covariates 
    vector[N_cultures] temp_c;
    vector[N_cultures] nutrient_c;

    // Outcome 
    vector[N] bacterial_count;    
}
// Data Transformation 
transformed data {
    // Center the covarites 
    vector[N_cultures] c_temp_c;
    vector[N_cultures] c_nutrient_c;
    c_temp_c = temp_c - mean(temp_c);
    c_nutrient_c = nutrient_c - mean(nutrient_c);
}
// Paramters 
parameters{
    // Paramters for the ODE intial state
    real alpha_0;
    real <lower = 0.001> omega_0;
    vector[N_cultures] z_0;

    // Paramters for the ODE rate term
    real alpha_k;
    real beta_temp;
    real beta_nutrient_c;
    real <lower = 0.001> omega_k;
    vector[N_cultures] z_k;

    // Observation Model paramter
    real <lower = 0.001> sigma_obs;
}
// Parameter Transformation 
transformed parameters {
    // Construct the ODE paramters as with linear covariate predictors 
    vector[N_cultures] log_co; // Initial State
    vector[N_cultures] log_r;  // Rate
    vector[N_cultures] co; 
    vector[N_cultures] r; 
    vector[N] log_mu;
    for(i in 1:N_cultures){
        log_co[i] = alpha_0 + omega_0 * z_0[i];
        log_r[i]  = alpha_k + beta_temp * c_temp_c[i] + beta_nutrient_c * c_nutrient_c[i] + omega_k * z_k[i];
    }
    co = exp(log_co);
    r  = exp(log_r);

    // Linear Predictor 
    for (i in 1:N)
        log_mu[i] = log_co[cultures_id[i]] + r[cultures_id[i]] * day_id[i];
}
// Model Block 
model{
    // Priors for the ODE intial state
    alpha_0 ~ normal(log(50), 1);
    omega_0 ~ exponential(1);
    z_0 ~ std_normal();

    // Priors for the ODE rate term
    alpha_k ~ normal(log(0.15), 1);
    beta_temp ~ normal(0.03, 0.1);
    beta_nutrient_c ~ normal(0.15, 0.01);
    omega_k ~ exponential(1);
    z_k ~ std_normal();

    // Priors observation Model paramter
    sigma_obs ~ exponential(1);

    // Model likelihood 
    if(prior_only == 0){
        // Obsevation Model 
        bacterial_count ~ lognormal(log_mu, sigma_obs);
    }

}
// Minimal Generated Quantities 
generated quantities {
    vector[N] log_lik;
    vector[N] y_rep;
    for (i in 1:N) {
        log_lik[i] = lognormal_lpdf(bacterial_count[i] | log_mu[i], sigma_obs);
        y_rep[i]   = lognormal_rng(log_mu[i], sigma_obs);
    }
}