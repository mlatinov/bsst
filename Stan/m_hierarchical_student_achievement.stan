// Input Data Layer 
data{
    // Outcome and prior switch 
    int<lower = 1> N;
    int<lower = 0, upper =1> prior_only;
    vector[N] test_scores;

    // Covariates 
    vector[N] ses_index;
    vector[N] class_size;
    vector[N] study_h;

    // Indexes 
    array[N] int school_id;
    array[N] int district_id;
    array[N] int state_id;
    int n_school_j;
    int n_district_k;
    int n_state_l;
}
// Model Paramters 
parameters{
    // Intercepts 
    real alpha;
    real <lower = 0.001> sigma_school;
    real <lower = 0.001> sigma_district;
    real <lower = 0.001> sigma_state;
    vector[n_school_j]   z_alpha_school;
    vector[n_district_k] k_alpha_district;
    vector[n_state_l]    l_alpha_state;
    
    // Multilevel coefficients 
    real ses_index_bar;
    real <lower = 0.001> sigma_ses_index;
    vector[n_state_l] ses_zl;
    
    // Fixed coefs 
    real beta_study_h;
    real beta_class_size;

    // Test Score sigma
    real <lower = 0.001> test_score_sigma;
}
// Parameter Transformations 
transformed parameters {
    // Transform the intercepts and the SES random Effect Non Centered 
    vector[n_school_j]   alpha_school;
    vector[n_district_k] alpha_district;
    vector[n_state_l]    alpha_state;
    vector[n_state_l]    beta_ses_l;
    for(i in 1:n_school_j){
        alpha_school[i] = sigma_school * z_alpha_school[i];
    }
    for(i in 1:n_district_k){
        alpha_district[i] = sigma_district * k_alpha_district[i];
    }
    for(i in 1:n_state_l){
        alpha_state[i] = sigma_state * l_alpha_state[i];
        beta_ses_l[i]   = ses_index_bar + sigma_ses_index * ses_zl[i];
    }
}
// Model Layer 
model{
    // Priors
    alpha         ~ normal(70, 2);
    ses_index_bar ~ normal(6, 1);
    sigma_district ~ exponential(3);
    sigma_school   ~ exponential(5);
    sigma_state    ~ exponential(2); 
    sigma_ses_index ~ exponential(1);
    test_score_sigma ~ exponential(1);
    z_alpha_school  ~ std_normal();
    k_alpha_district ~ std_normal();
    l_alpha_state    ~ std_normal();
    ses_zl           ~ std_normal();
    beta_study_h    ~ normal(4, 1);
    beta_class_size ~ normal(-0.3, 1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            // Calculate Mu 
            mu[i] =
                alpha 
                + alpha_state[state_id[i]]
                + alpha_district[district_id[i]]
                + alpha_school[school_id[i]]
                + beta_ses_l[state_id[i]] * ses_index[i]
                + beta_study_h * study_h[i]
                + beta_class_size * class_size[i];
        }
        // Sample Test Scores form Normal Distribition 
        test_scores ~ normal(mu, test_score_sigma);
    }
}
// Minimal generated quantites 
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i] =
            alpha 
            + alpha_state[state_id[i]]
            + alpha_district[district_id[i]]
            + alpha_school[school_id[i]]
            + beta_ses_l[state_id[i]] * ses_index[i]
            + beta_study_h * study_h[i]
            + beta_class_size * class_size[i];

        log_lik[i] =
            normal_lpdf(
                test_scores[i] | mu[i], test_score_sigma
            );
        y_rep[i] =
            normal_rng(mu[i], test_score_sigma);
    }
}