// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;
    int <lower = 2> N_installations;
    array[N] int  installation_id;

    // Spline Settings 
    int <lower = 2> N_basis;
    matrix[N, N_basis] B_temperature;

    // Covariates 
    vector[N] dust_index;

    // Outcome 
    vector[N] efficiency;
}

// Model Parameters
parameters{
    // Parameters for Correlated Random Intercept and Beta Dust index
    real alpha_bar;
    real beta_dust_bar;
    vector <lower = 0.001>[2] sds;
    cholesky_factor_corr[2] L;
    matrix[2, N_installations] z_matrix;

    // Parameters for the P-Spline over the Temperature
    real wk_pop1;
    real wk_pop2;
    real <lower = 0.001> sd_wk_pop;
    vector[N_basis - 2] zk_pop;

    // NCP Parameters for the P-spline over the Temperature 
    matrix[N_basis, N_installations] z_unit;
    real <lower = 0.001> omega;

    // Observational Model Parameters 
    real <lower = 0.001> kappa;
}
// Transform Parameters
transformed parameters {
    // Ensemble Correlated Effects for the intercept and beta dust index coef
    matrix[2, N_installations] v;
    v = diag_pre_multiply(sds, L) * z_matrix;

    // Ensemble weights Temperature P-Spline 
    vector[N_basis] wk;
    wk[1] = wk_pop1;
    wk[2] = wk_pop2;
    for(k in 3:N_basis){
        wk[k] = 2 * wk[k - 1] - wk[k - 2] + sd_wk_pop * zk_pop[k - 2];
    }
    
    // P-spline NCP 
    matrix[N_basis, N_installations] w_units;
    for(i in 1:N_installations){
        w_units[, i] = wk + omega .* z_unit[, i];
    } 
}
// Model 
model{
    // Priors for Correlated Random Intercept and Beta Dust index
    alpha_bar ~ normal(0, 1);
    beta_dust_bar ~ normal(0, 1);
    sds ~ exponential(1);
    L   ~ lkj_corr_cholesky(2);
    to_vector(z_matrix) ~ std_normal();

    // Priors for the P-Spline over the Temperature
    wk_pop1 ~ normal(0, 1);
    wk_pop2 ~ normal(0, 1);
    sd_wk_pop ~ exponential(1);
    zk_pop    ~ std_normal();

    // NCP Priors for the P-spline over the Temperature 
    to_vector(z_unit) ~ std_normal();
    omega ~ exponential(1);

    // Observational Model Priors 
    kappa ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = inv_logit( 
                (alpha_bar + v[1, installation_id[i]]) 
                +
                (beta_dust_bar + v[2, installation_id[i]]) * dust_index[i]
                +
                (B_temperature[i] * w_units[, installation_id[i]])
            );
        }
    // Observation Beta Model Model 
    efficiency ~ beta(mu * kappa, (1 - mu) *  kappa);
    }

}
// Minimal Generated Quantities
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i] = inv_logit(
            (alpha_bar + v[1, installation_id[i]]) 
            +
            (beta_dust_bar + v[2, installation_id[i]]) * dust_index[i]
            +
            (B_temperature[i] * w_units[, installation_id[i]])
        ); 
        log_lik[i] = beta_lpdf(efficiency[i] | mu[i] * kappa, (1 - mu[i]) * kappa);
        y_rep[i]   = beta_rng(mu[i] * kappa, (1 - mu[i]) * kappa);
    }
   
}