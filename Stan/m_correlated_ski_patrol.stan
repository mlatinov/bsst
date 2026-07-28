// Input data layer 
data{
    // Number of observartions and prior only switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] briefings_h;

    // Indexes
    int <lower = 2> n_resorts;
    array[N] int resort_id;

    // Outcome 
    array[N] int incident;
}

// Model Paramters 
parameters{
    // Correlated intercept and briefing slope 
    // Population intecept and population briefing slope 
    real a;
    real b;
    
    // Sd of a and b 
    vector <lower = 0.001>[2] sigmas;
    
    // Cholesky Factor 
    cholesky_factor_corr[2] L;
    matrix[2, n_resorts] z_matrix;
}

// Transform the correlation paramters
transformed parameters {
   matrix[2, n_resorts] v;
   v = diag_pre_multiply(sigmas, L) * z_matrix;
}
// Model 
model{
    // Correlated intercept and briefing slope priors 
    a ~ normal(0.003, 0.001);
    b ~ normal(-2, 1); 
    sigmas ~ exponential(1);
    L      ~ lkj_corr_cholesky(2) ;
    to_vector(z_matrix) ~ std_normal();

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] pi;
        for(i in 1:N){
            pi[i] = inv_logit(
                (a + v[1, resort_id[i]]) + (b + v[2, resort_id[i]]) * briefings_h[i]
            );
        }
    // Observational Model Sampled from Bernoulli Distribution 
    incident ~ bernoulli(pi);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] pi;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        pi[i] = inv_logit(
                (a + v[1, resort_id[i]]) + (b + v[2, resort_id[i]]) * briefings_h[i]
        );
        log_lik[i] = bernoulli_lpmf(incident[i] | pi[i]);
        y_rep[i]   = bernoulli_rng(pi[i]);
    }
    // Correlation Recovery Specific for the Correlated Models 
    matrix[2,2] Rho = multiply_lower_tri_self_transpose(L); 
}