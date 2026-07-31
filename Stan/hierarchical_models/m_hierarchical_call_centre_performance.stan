// Data input layer 
data{
    // Num of Observations and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] call_min;
    vector[N] complexity;
    vector[N] agent_year_experience;
    vector[N] repeated_call;

    // Indexes 
    int <lower = 2> n_call_centers;
    array[N] int call_center_id;
    
    // Outcome 
    array[N] int<lower=0, upper=1> resolved_problem;
}
// Model Paramters 
parameters{
    // Random Effects 
    real baseline_resolution;
    real sigma_resolution;
    vector[n_call_centers] zj;
    
    // Fixed effects 
    real beta_call_min;
    real beta_compexity;
    real beta_agent_year_experience;
    real beta_repeated_call;
}
// Random Effect Parameter Transformation Non-Centered
transformed parameters {
   vector[n_call_centers] alpha_j;
   for(i in 1:n_call_centers){
    alpha_j[i] = baseline_resolution + sigma_resolution * zj[i];
   }
}
// Model 
model{
    // Random Effects Priors 
    baseline_resolution ~ normal(1.2, 0.1);
    sigma_resolution ~ exponential(1);
    zj ~ std_normal();

    // Fixed effects Priors 
    beta_call_min  ~ normal(-0.08, 0.01);
    beta_compexity ~ normal(-0.60, 0.1);
    beta_agent_year_experience ~ normal(0.30, 0.1);
    beta_repeated_call         ~ normal(-1.10, 1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] pi;
        for(i in 1:N){
            pi[i] = inv_logit(
                alpha_j[call_center_id[i]]
                + beta_call_min  * call_min[i]
                + beta_compexity * complexity[i]
                + beta_repeated_call * repeated_call[i]
                + beta_agent_year_experience *  agent_year_experience[i]
            );
        }
    // Observation Model 
    resolved_problem ~ bernoulli(pi);   
    }
}
// Minimal Generated Quantites
generated quantities {
    vector[N] pi;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
         pi[i] = inv_logit(
                alpha_j[call_center_id[i]]
                + beta_call_min  * call_min[i]
                + beta_compexity * complexity[i]
                + beta_repeated_call * repeated_call[i]
                + beta_agent_year_experience *  agent_year_experience[i]
            );
        log_lik[i] = bernoulli_lpmf(resolved_problem[i] | pi[i]);
        y_rep[i]   = bernoulli_rng(pi[i]);
    }
}