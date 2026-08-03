// Data Input Layer 
data{
    // Number of observations and prior only switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Indexes
    int n_bakery;
    array[N] int bakery_id;

    // Covariates 
    vector[N] temperature_c;
    vector[N] is_weekend;
    vector[N] foot_trafic;

    // Outcome 
    vector[N] daily_revenue;
}
// Transform Data 
transformed data {
   // Standartize numerical covariates 
   vector[N] z_temperature;
   vector[N] z_foot_trafic;
   z_temperature = (temperature_c - mean(temperature_c)) / sd(temperature_c);
   z_foot_trafic = (foot_trafic - mean(foot_trafic)) / sd(foot_trafic);
}
// Model Paramters 
parameters{
    // Random Intercept Paramters 
    real z_baseline_daily_revenue;
    real <lower = 0.001> bakery_revenue_sd;
    vector[n_bakery] zj;

    // Fixed Effects Paramters 
    real z_beta_temp;
    real z_beta_trafic;
    real beta_weekend;

    // Outcome paramters 
    real <lower = 0.001> sd_revenue; 
}
// Transform Paramters for the Random Intercept
transformed parameters {
   vector[n_bakery] alpha_j;
   for(i in 1:n_bakery){
    alpha_j[i] = z_baseline_daily_revenue + bakery_revenue_sd * zj[i];
   }
}
// Model 
model{
    // Random Intercept paramters 
    z_baseline_daily_revenue ~ normal(100, 5);
    bakery_revenue_sd ~ exponential(1);
    zj ~ std_normal();

    // Fixed Effects Priors 
    z_beta_temp    ~ normal(0, 1);
    z_beta_trafic  ~ normal(0, 1);
    beta_weekend ~ normal(0, 1);

    // Outcome sd 
    sd_revenue ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] mu;
        for(i in 1:N){
            mu[i] = 
                alpha_j[bakery_id[i]] 
                + z_beta_temp * z_temperature[i]
                + z_beta_trafic * z_foot_trafic[i]
                + beta_weekend * is_weekend[i]; 
        }
    // Outcome model 
    daily_revenue ~ normal(mu, sd_revenue);
    }
}
// Minimal Generated Quantites 
generated quantities { 
    // natural-scale quantities
    real beta_temp    = z_beta_temp  / sd(temperature_c);
    real beta_trafic  = z_beta_trafic / sd(foot_trafic);
    real baseline_daily_revenue = z_baseline_daily_revenue - beta_temp  * mean(temperature_c) - beta_trafic  * mean(foot_trafic) ;

    // Reproduce the variables 
    vector[N] mu;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        mu[i] = 
            alpha_j[bakery_id[i]] 
            + z_beta_temp * z_temperature[i]
            + z_beta_trafic * z_foot_trafic[i]
            + beta_weekend * is_weekend[i]; 
        log_lik[i] = normal_lpdf(daily_revenue[i] | mu[i], sd_revenue);
        y_rep[i]   = normal_rng(mu[i], sd_revenue);
    }
}