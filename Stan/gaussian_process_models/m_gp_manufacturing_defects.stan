// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 2> N_factory;
    int <lower = 0, upper = 1> prior_only;
    array[N] int factory_id;

    // Covariates 
    vector[N] speed;

    // Outcome 
    array[N] int defects_count;
}
// Transform Data
transformed data {
    // Small Delta constant for the cholessky decomposition 
    real delta = 1e-9;

    // Count how many observations belogns to each factory
    array[N] int factory_size = rep_array(0, N_factory);
    for(i in 1:N){
        factory_size[factory_id[i]] += 1;
    }

}
// Model Parameters
parameters{
    // GP Hierarchical paramters for rho controling the lenght scale
    real mean_rho;
    real <lower = 0.001> sd_rho;
    vector[N_factory] z_rho;

    // GP Hierarchical parameters for alpha controling the amplitude 
    real mean_alpha;
    real <lower = 0.001> sd_alpha;
    vector[N_factory] z_alpha;

    // Raw Gaussian vector 
    vector[N] z;

}
// Transform Parameters
transformed parameters {
    // Recover GP Hierarchical parameters 
    vector[N_factory] rho;
    vector[N_factory] alpha;
    rho = exp(mean_rho + sd_rho * z_rho);
    alpha = exp(mean_alpha + sd_alpha * z_alpha);

    // Latent GP function of speed 
    vector[N] f_speed;

    // Process every factory seperatly 
    for(j in 1:N_factory){
        // Count the number of observations in each 
        int nj = factory_size[j];

        // Store speed values  
        array[nj] real x_j;
        vector[nj] zj;
        int idx = 1;

        // Extract the observations belonging to j factory 
        for(i in 1:N){
            if(factory_id[i] == j){
                x_j[idx] = speed[i];
                zj[idx]  =  z[i];
                idx += 1;
            }
        }

        // Build the covariance matrix K
        matrix[nj, nj] K = add_diag(gp_exp_quad_cov(x_j, alpha[j], rho[j]), delta);
        
        // Build the Transformation matrix from the covariance matrix K
        matrix[nj, nj] L = cholesky_decompose(K);

        // Create a smooth GP curve 
        vector[nj] f =  L * z;

        // Return the GP values back in function of speed vector 
        idx = 1;
        for(i in 1:N){
            if(factory_id[i] == j){
                f_speed[i] = f[idx];
                idx += 1;
            }
        }
    }
}

// Model 
model{
    // GP Hierarchical priors for rho controling the lenght scale
    mean_rho ~ normal(0, 1);
    sd_rho ~ exponential(1);
    z_rho  ~ std_normal();

    // GP Hierarchical priors for alpha controling the amplitude 
    mean_alpha ~ normal(0, 1);
    sd_alpha   ~ exponential(1);
    z_alpha    ~ std_normal();

    // Raw Gaussian vector priors
    z ~ std_normal();

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] lambda;
        for(i in 1:N){
            lambda[i] = exp(f_speed[i]);
        }
    // Observation model 
    defects_count ~ poisson(lambda);
    }
}
// Minimal Generated Quantities
generated quantities {
    vector[N] lambda;
    vector[N] defects_count_rep;
    vector[N] log_lik;
    for(i in 1:N){
        lambda[i] = exp(f_speed[i]);
        defects_count_rep[i] = poisson_lpmf(defects_count[i] | lambda[i]);
        log_lik[i] = poisson_rng(lambda[i]);
    }
}