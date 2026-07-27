// Data input Layer 
data{
    // N and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;
    
    // Indexes
    array[N] int shop_id;
    int <lower = 2> n_shops;
    
    // Covariates 
    vector[N] experiance;
    vector[N] training_h;

    // Outcome 
    vector[N] rating;
}
// Model Paramters 
parameters{
    // Parameters for varing shop id intercept 
    real alpha_j_bar;
    real <lower = 0.001> shop_sigma;
    vector[n_shops] zj;

    // Paramters for fixed covariates effects 
    real beta_training_h;
    real beta_experiance;

    // Outcome Sd 
    real <lower = 0.001> sd_rating; 
}
// Paramter transformation for the varing intercept 
transformed parameters {
   vector[n_shops] alpha_j;
   for(i in 1:n_shops){
    alpha_j[i] = alpha_j_bar + shop_sigma * zj[i]; 
   }
}
// Model 
model{
    // Priors 
    alpha_j_bar ~ normal(5.5, 0.5);
    shop_sigma  ~ exponential(0.5);
    sd_rating   ~ exponential(1);
    zj          ~ std_normal();
    beta_training_h ~ normal(0.1, 0.1);
    beta_experiance ~ normal(0.2, 1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = 
                alpha_j[shop_id[i]] 
                + beta_experiance * experiance[i] 
                + beta_training_h * training_h[i];
        }
        rating ~ normal(mu, sd_rating);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
         mu[i] = 
            alpha_j[shop_id[i]] 
            + beta_experiance * experiance[i] 
            + beta_training_h * training_h[i];
        log_lik[i] =
            normal_lpdf(
                rating[i] | mu[i], sd_rating
            );
        y_rep[i] =
            normal_rng(mu[i], sd_rating);
    }
}