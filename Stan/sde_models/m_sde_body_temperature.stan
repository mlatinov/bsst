// Data input Layer 
data{
    // N Observations and Indexing 
    int <lower = 2> dt;
    int <lower = 2> N_hours;
    int <lower = 2> N_subjects;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N_subjects] ambient_temp;
    vector[N_subjects] activity_level;
    
    // Outcome 
    matrix[N_hours, N_subjects] body_temp;
}
// Transform Data 
transformed data {
    // Center the covariates 
    vector[N_subjects] c_ambient_temp;
    vector[N_subjects] c_activity_level;
    c_ambient_temp   = ambient_temp - mean(ambient_temp);
    c_activity_level = activity_level - mean(activity_level); 
}

// Model Paramters 
parameters{
    // Heat inflow paramters for a
    real mean_heat_inflow_rate;
    real beta_ambient_t;
    real <lower = 0.001> sd_subjects_heat_inflow_rate;
    vector[N_subjects] z_a;

    // Heat outflow parameters for b 
    real mean_outflow_rate;
    real beta_activity;
    real <lower = 0.001> sd_subjects_outflow_rate;
    vector[N_subjects] z_b;

    // Initial State Paramter
    real mean_intial_state;
    real <lower = 0.001> sd_initial_state;
    vector[N_subjects] z_in;

    // Paramteters for the OU process
    real <lower = 0.001> sd_ou;
    matrix[N_hours, N_subjects] C_raw;

    // Observation parameters 
    real <lower = 0.001> sd_obs;
}
// Transform Model Paramteters 
transformed parameters {
    // Ensemble the paramters a and b and intiial state xo  
    vector[N_subjects] a;
    vector[N_subjects] b;
    vector[N_subjects] xo;
    for(i in 1:N_subjects){
        a[i] = exp(mean_heat_inflow_rate + beta_ambient_t * c_ambient_temp[i] + sd_subjects_heat_inflow_rate * z_a[i]);
        b[i] = exp(mean_outflow_rate + beta_activity * c_activity_level[i] + sd_subjects_outflow_rate * z_b[i]);
        xo[i] = mean_intial_state + sd_initial_state * z_in[i];
    }

    // OU Process 
    matrix[N_hours, N_subjects] C;
    
    // Initiliaze the process 
    C[1, ] = to_row_vector(xo);
    
    // Continue the OU Process 
    for(t in 2:N_hours){
        vector[N_subjects] mean_t = a ./ b + (to_vector(C[t - 1,]) - a ./ b) .* exp(-b * dt);
        vector[N_subjects] sd_t   = sqrt((sd_ou^2 / (2 * b)) .* (1 - exp(-2 * b * dt)));
        C[t, ] = to_row_vector(mean_t + sd_t .* to_vector(C_raw[t, ]));
    } 
   
}
// Model 
model{
    // Heat inflow priors for a
    mean_heat_inflow_rate ~ normal(log(0.15) + log(37), 0.1);
    beta_ambient_t ~ normal(0.015, 0.01);
    sd_subjects_heat_inflow_rate ~ exponential(1);
    z_a ~ std_normal();

    // Heat outflow priors for b 
    mean_outflow_rate ~ normal(log(0.15), 0.1);
    beta_activity ~ normal(-0.03, 0.01);
    sd_subjects_outflow_rate ~ exponential(1);
    z_b ~ std_normal();

    // Initial State priors
    mean_intial_state ~ normal(36.8, 0.5);
    sd_initial_state  ~ exponential(1);
    z_in ~ std_normal();

    // Priors for the OU process
    sd_ou ~ exponential(1);
    to_vector(C_raw) ~ std_normal();

    // Observation priors 
    sd_obs ~ exponential(1);

    // Observational Model 
    if(prior_only  == 0){
        to_vector(body_temp) ~ normal(to_vector(C), sd_obs);
    }
}
// Minimal Generated quantities 
generated quantities{
    matrix[N_hours, N_subjects] body_temp_rep;
    matrix[N_hours, N_subjects] log_lik;
    for(i in 1:N_hours){
        for(j in 1:N_subjects){
            body_temp_rep[i, j] = normal_rng(C[i, j], sd_obs);
            log_lik[i, j] = normal_lpdf(body_temp[i, j] | C[i, j], sd_obs);
        }
    }
}