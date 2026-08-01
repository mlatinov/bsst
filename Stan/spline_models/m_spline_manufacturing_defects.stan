// Data Input Layer 
data{
    // Settings 
    int <lower = 1> N;
    int <lower = 2> N_factory;
    int <lower = 0, upper = 1> prior_only;
    array[N] int factory_id;

    // Spline Settings 
    int <lower = 2> N_basis;
    matrix[N, N_basis] B_speed;

    // Outcome 
    array[N] int defects_count;  
}
// Model Parameters
parameters{
    // Random Intercept Paramters 
    real mean_factories_defects;
    real <lower = 0.001> sd_factories_defects;
    vector[N_factory] z_a;

    // P-Splines Population Parameters
    real wk_pop1;
    real wk_pop2;
    real <lower = 0.001> sd_wk_pop;
    vector[N_basis - 2] z_wk_pop;

    // Per Unit P-Splines Parameters
    matrix[N_basis, N_factory] z_unit;
    real <lower = 0.001> omega;
}

// Transform Parameters
transformed parameters {
   // Ensemble Random intercept NCP
   vector[N_factory] alpha_j;
   for(i in 1:N_factory){
    alpha_j[i] = mean_factories_defects + sd_factories_defects * z_a[i];
   }

   // Ensemble P-spline NCP
   vector[N_basis] wk_pop;
   wk_pop[1] = wk_pop1;
   wk_pop[2] = wk_pop2;
   for(k in 3:N_basis){
    wk_pop[k] = 2 * wk_pop[k - 1] - wk_pop[k - 2] + sd_wk_pop * z_wk_pop[k - 2];
   }
   
   // NCP
   matrix[N_basis, N_factory] w_unit;
   for(i in 1:N_factory){
    w_unit[,i] = wk_pop + omega  .* z_unit[,i] ;
   }
}
// Model 
model{
    // Random Intercept Priors 
    mean_factories_defects ~ normal(0, 1);
    sd_factories_defects ~ exponential(1);
    z_a ~ std_normal();

    // P-Splines Population Priors
    wk_pop1 ~ normal(0, 1);
    wk_pop2 ~ normal(0, 1);
    sd_wk_pop ~ exponential(1);
    z_wk_pop  ~ std_normal();

    // Per Unit P-Splines Priors
    to_vector(z_unit) ~ std_normal() ;
    omega ~ exponential(1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] lambda;
        for(i in 1:N){
            lambda[i] = exp(alpha_j[factory_id[i]] + B_speed[i] * z_unit[, factory_id[i]]);
        }
    // Observation model 
    defects_count ~ poisson(lambda);
    }
}
// Minimal Generated Quantities
generated quantities {
    vector[N] lambda;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        lambda[i] = alpha_j[factory_id[i]] + B_speed[i] * w_unit[, factory_id[i]];
        log_lik[i] = poisson_lpmf(defects_count[i] | lambda[i]);
        y_rep[i]   = poisson_rng(lambda[i]);
    }
   
}