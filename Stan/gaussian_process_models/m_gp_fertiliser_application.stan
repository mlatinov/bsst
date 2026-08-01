// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 2> N_soil_type;
    int <lower = 0, upper = 1> prior_only;
    array[N] int soil_type_id;

    // Covariates 
    vector[N] fertiliser;

    // Outcome 
    vector[N] biomass;
}
// Transform Data
transformed data {
    // Delta ...Cholesky decomposition succeeds even when K is nearly singular.
    real delta = 1e-9;

    // Count how many observations belong to each soil type
    array[N_soil_type] int soil_type_size = rep_array(0, N_soil_type);
    for(i in 1:N){
        soil_type_size[soil_type_id[i]] += 1;
    }

}
// Model Parameters
parameters{
    // Hierarchical parameters for GP length-scales (rho)
    real mean_rho;
    real <lower = 0.001> tau_rho;
    vector[N_soil_type]  z_rho;

    // Hierarchical parameters for GP amplitudes (alpha)
    real mean_alpha;
    real <lower = 0.001> tau_alpha;
    vector[N_soil_type]  z_alpha;

    // Raw latent Gaussian vectors parameters
    vector[N] z;

    // Observational Model Parameters
    real <lower = 0.001> sd_obs;
}

// Transform Parameters
transformed parameters {
   // Recover each unit's GP hyperparameters
   vector[N_soil_type] rho;
   vector[N_soil_type] alpha;
   
   rho = exp(mean_rho + tau_rho * z_rho);
   alpha = exp(mean_alpha + tau_alpha * z_alpha);

   // Latent GP functions
   vector[N] f_fertilizer;

   // Process every soil type separately.
   for(j in 1:N_soil_type){
    
   // Number of observations for this soil type
    int nj = soil_type_size[j];

    // Store fertilizer values for this soil only and raw values 
    array[nj] real x_j;
    vector[nj] z_j;
    int idx = 1;

    // Extract observations belonging to soil j
    for(i in 1:N){
        if(soil_type_id[i] == j){ 
            x_j[idx] = fertiliser[i];
            z_j[idx] = z[i];
            idx += 1;       
        }
    }

    // Build covariance matrix
    matrix[nj,nj] K_j = add_diag(gp_exp_quad_cov(x_j, alpha[j], rho[j]), delta);

    // Convert covariance matrix into a transformation matrix
    matrix[nj,nj] L_j = cholesky_decompose(K_j);

    // Create smooth GP curve
    vector[nj] f_j = L_j * z_j;
    
    // Put the GP values back into their original rows
    idx = 1;
    for(i in 1:N){
        if(soil_type_id[i] == j){
            f_fertilizer[i] = f_j[idx];
            idx += 1;
        }
    }
   }
}

// Model 
model{
    // Hierarchical prior for GP length-scales (rho)
    mean_rho ~ normal(0, 1);
    tau_rho  ~ exponential(1);
    z_rho    ~ std_normal();

    // Hierarchical prior for GP amplitudes (alpha)
    mean_alpha ~ normal(0, 1);
    tau_alpha  ~ exponential(1);
    z_alpha    ~ std_normal();

    // Raw latent Gaussian vectors prior
    to_vector(z) ~ std_normal();

    // Observational Model prior
    sd_obs ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        // Observation model 
        biomass ~ normal(f_fertilizer, sd_obs);
    }
}
// Minimal Generated Quantities
generated quantities {
    vector[N] biomass_rep;
    vector[N] log_lik;
    for(i in 1:N){
        biomass_rep[i] = normal_rng(f_fertilizer[i], sd_obs);
        log_lik[i] = normal_lpdf(biomass[i] | f_fertilizer[i], sd_obs);
    }
}