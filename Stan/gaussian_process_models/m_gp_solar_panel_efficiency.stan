// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 2> N_installations;
    int <lower = 0, upper = 1> prior_only;
    array[N] int installation_id;

    // Covariates 
    vector[N] temperature;
    vector[N] dust_index;

    // Outcome 
    vector[N] efficiency;
}
// Transform Data
transformed data {
    // Center the temperature 
    vector[N] c_temperature; 
    c_temperature = temperature - mean(temperature);

    // Count the observation in each installation 
    array[N] int installation_size = rep_array(0, N_installations);
    for(i in 1:N_installations){
        installation_size[installation_id[i]] += 1;
    }  
}
// Model Parameters
parameters{
    // GP Hierarchical parameters for rho controling the lenght scale 
    real mean_gp_rho;
    real <lower = 0.001> sd_gp_rho;
    vector[N_installations] z_gp_rho;

    // GP Hierarchical parameters for alpha controling the amplitude
    real mean_gp_alpha;
    real <lower = 0.001> sd_gp_alpha;
    vector[N_installations] z_gp_alpha;

    // Raw Gaussian Vector
    vector[N] z_gp;

    // Correlated Hierarchical intercept and beta dust parameters 
    real alpha_bar;
    real beta_dust_var;
    vector <lower = 0.001> [2] sds;
    cholesky_factor_corr[2] L;
    matrix[2, N_installations] z_cor;

    // Observational Model Parameters  
    real <lower = 0.001> kappa;

}
// Transform Parameters
transformed parameters {
    // Recover the Correlated Hierarchical Parameters for the intercept and beta dust index 
    matrix[2, N_installations] v;
    v = diag_pre_multiply(sds, L) * z_cor;

    // Recover the GP rho and alpha parameters 
    vector[N_installations] gp_rho;
    vector[N_installations] gp_alpha;
    gp_rho = exp(mean_gp_rho + sd_gp_rho * z_gp_rho);
    gp_alpha = exp(mean_gp_alpha + sd_gp_alpha * z_gp_alpha); 

    // Small Stability delta for the K convariace matrix 
    real delta = 1e-9;

    // Latent GP functon of temperature 
    vector[N] f_temperature;

    // Process every installation seperatly 
    for(j in 1:N_installations){

        // Count the number of observations in each
        int nj = installation_size[j];

        // Store speed values  
        array[nj] real x_j;
        vector[nj] zj;
        int idx = 1;

        // Extract the observations belonging to j installation 
        for(i in 1:N){
            if(installation_id[i] == j){
                x_j[idx] = c_temperature[i];
                zj[idx]  =  z_gp[i];
                idx += 1;
            }
        }

        // Build and transform the covariance matrix K 
        matrix[nj, nj] K = add_diag(gp_exp_quad_cov(x_j, gp_alpha[j], gp_rho[j]), delta);
        matrix[nj, nj] L_k = cholesky_decompose(K);

        // Create a smooth GP curve 
        vector[nj] f = L_k * z_gp;

        // Return the GP values back in function of speed vector 
        idx = 1;
        for(i in 1:N){
            if(installation_id[i] == j){
                f_temperature[i] = f[idx];
                idx += 1;
            }
        }
    }
}

// Model 
model{
    // GP Hierarchical priors for rho controling the lenght scale 
    mean_gp_rho ~ normal(0, 1);
    sd_gp_rho   ~ exponential(1) ;
    z_gp_rho ~ std_normal();

    // GP Hierarchical priors for alpha controling the amplitude
    mean_gp_alpha ~ normal(0, 1);
    sd_gp_alpha ~ exponential(1);
    z_gp_alpha ~ std_normal();

    // Raw Gaussian Vector
    z_gp ~ std_normal();

    // Correlated Hierarchical intercept and beta dust priors 
    alpha_bar ~ normal(0, 1);
    beta_dust_var ~ normal(0, 1);
    sds           ~ exponential(1);
    L             ~ lkj_corr_cholesky(2);
    to_vector(z_cor) ~ std_normal();

    // Observational Model priors  
    kappa ~ exponential(1);

    // Model Likelihood 
    if(prior_only  == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = inv_logit(
                (alpha_bar + v[1, installation_id[i]])                     // Random Intercept
                +
                (beta_dust_var + v[2, installation_id[i]]) * dust_index[i] // Random Dust Slope 
                +
                f_temperature[i]                                           // GP temperature function 
            );
        }
    // Observational Model 
    efficiency ~ beta(mu * kappa, (1 - mu) * kappa);
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
            (beta_dust_var + v[2, installation_id[i]]) * dust_index[i]
            +
            f_temperature[i]
        ); 
        log_lik[i] = beta_lpdf(efficiency[i] | mu[i] * kappa, (1 - mu[i]) * kappa);
        y_rep[i]   = beta_rng(mu[i] * kappa, (1 - mu[i]) * kappa);
    }
}