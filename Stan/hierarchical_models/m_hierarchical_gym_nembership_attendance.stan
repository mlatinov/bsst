// Data Input Layer 
data{
    // Number of observation and prior only switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] age;
    vector[N] distance;
    vector[N] trainer;

    // Indexes
    int <lower = 2> n_gyms;
    array [N] int gym_id;

    // Outcome 
    array[N] int visits;
}
// Model Paramters 
parameters{
    // Random Intercet
    real mean_log_baseline;
    real <lower = 0.001> sigma_gym_visits;
    vector[n_gyms] zj_a;
    
    // Random age effect 
    real mean_age_effect;
    real <lower = 0.001> sigma_age_effects;
    vector[n_gyms] zj_age;

    // Random Distance Effect 
    real mean_distance_effect;
    real <lower = 0.001> sigma_distance_effects;
    vector[n_gyms] zj_d;

    // Random Trainer Effect 
    real mean_trainer_effect;
    real <lower = 0.001> sigma_trainer_effects;
    vector[N] zj_t;
}
// Random Effects Transformations Paramters
transformed parameters {
   vector[n_gyms] alpha_j;
   vector[n_gyms] beta_age_j;
   vector[n_gyms] beta_distance_j;
   vector[n_gyms] beta_trainer_j;
   for(i in 1:n_gyms){
    alpha_j[i]    = mean_log_baseline + sigma_gym_visits * zj_a[i];
    beta_age_j[i] = mean_age_effect   + sigma_age_effects * zj_age[i];
    beta_distance_j[i] = mean_distance_effect + sigma_distance_effects * zj_d[i];
    beta_trainer_j[i]  = mean_trainer_effect  + sigma_trainer_effects  * zj_t[i];
   }
}
// Model 
model{
    // Random Intercet Priors
    mean_log_baseline ~ normal(log(10), 1);
    sigma_gym_visits ~ exponential(1);
    zj_a ~ std_normal();
    
    // Random age effect Priors 
    mean_age_effect ~ normal(-0.01, 0.001);
    sigma_age_effects ~ exponential(1);
    zj_age ~ std_normal();

    // Random Distance Effect Priors 
    mean_distance_effect ~ normal(0.02, 0.01);
    sigma_distance_effects ~ exponential(1);
    zj_d ~ std_normal();

    // Random Trainer Effect Priors 
    mean_trainer_effect ~ normal(0.25, 0.02);
    sigma_trainer_effects ~ exponential(1);
    zj_t ~ std_normal();

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] lambda;
        for(i in 1:N){
            lambda[i] = exp(
                alpha_j[gym_id[i]] 
                + beta_age_j[gym_id[i]] * age[i]
                + beta_distance_j[gym_id[i]] * distance[i]
                + beta_trainer_j[gym_id[i]]  * trainer[i]
            );
        }
        // Observation Model 
        visits ~ poisson(lambda);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] lambda;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        lambda[i] = exp(
            alpha_j[gym_id[i]] 
            + beta_age_j[gym_id[i]] * age[i]
            + beta_distance_j[gym_id[i]] * distance[i]
            + beta_trainer_j[gym_id[i]]  * trainer[i]
        );
        log_lik[i] = poisson_lpmf(visits[i] | lambda[i]);
        y_rep[i]   = poisson_rng(lambda[i]);
    }
}
