// Data input Layer 
data{
    // Observation N and prior only switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] sun_hours;
    vector[N] panel_temp;
    vector[N] dust_index;

    // Indexes 
    int <lower = 2> n_farms;
    array[N] int farm_id;

    // Outcome 
    vector[N] kw;
}

// Model Paramters 
parameters{
    // Correlated Intercept, Sun hours slope, Panel Temp slope 
    real a; // Intercept
    real b; // Sun_h slope 
    real c; // Panel temp slope
    vector <lower = 0.001>[3]  sigmas;
    cholesky_factor_corr[3] L;
    matrix[3, n_farms] zj;

    // Dust index slope paramter and residual sigma
    real beta_dust;
    real <lower = 0.01> sigma;
}
// Transform correlated model paramters 
transformed parameters {
   matrix[3, n_farms] v;
   v = diag_pre_multiply(sigmas, L) * zj;
}
// Model 
model{
    // Correlated Intercept, Sun hours slope, Panel Temp slope Priors 
    a ~ normal(18, 1);      // Intercept
    b ~ normal(6, 1);       // Sun_h slope 
    c ~ normal(-0.75, 1);   // Panel temp slope
    sigmas        ~ exponential(1);
    L             ~ lkj_corr_cholesky(2);
    to_vector(zj) ~ std_normal();

    // Dust index slope paramter and residual sigma Priors
    beta_dust ~ normal(-12, 1);
    sigma     ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = 
                (a + v[1, farm_id[i]])                  // Intercept 
                + 
                (b + v[2, farm_id[i]]) * sun_hours[i]   // Sun Slope
                +
                (c + v[3, farm_id[i]]) * panel_temp[i]  // Panel Temp Slope
                + 
                beta_dust              * dust_index[i]; // Dust Slope 
        }
        // Obervation Model 
        kw ~ normal(mu, sigma);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i] = 
            (a + v[1, farm_id[i]])                  // Intercept 
            + 
            (b + v[2, farm_id[i]]) * sun_hours[i]   // Sun Slope
            +
            (c + v[3, farm_id[i]]) * panel_temp[i]  // Panel Temp Slope
            + 
            beta_dust              * dust_index[i];
        
        log_lik[i] = normal_lpdf(kw[i] | mu[i], sigma);
        y_rep[i]   = normal_rng(mu[i], sigma);
    }
    // Correlation Recovery Specific for the Correlated Models 
    matrix[3,3] Rho = multiply_lower_tri_self_transpose(L); 
}