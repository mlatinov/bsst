// Data input Layer 
data{
    // Observations and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] dose;
    vector[N] comorbidities;
    vector[N] age;

    // Indexing 
    int <lower = 2> N_hospitals;
    int <lower = 2> N_wards;
    array[N] int hospital_id;
    array[N] int ward_id;

    // Outcome 
    array[N] int recovery;
}
// Model Parameters 
parameters{
    // Paramters for Correlated hospital intercept and Dose effect 
    real a_bar;                         // Hospital intercept mean 
    real b_bar;                         // Dose Effect mean  
    vector <lower = 0.001 >[2] sigmas;
    matrix[2, N_hospitals]     z_matrix;
    cholesky_factor_corr[2]    L;

    // Ward intercept parameters
    real ward_a_bar;
    real <lower = 0.001> ward_sigma;
    vector[N_wards] zj;

    // Fixed Effects 
    real beta_comorbidities;
    real beta_age;
}
// Transform paramters 
transformed parameters {
    // Compose the correlated effects 
    matrix[2, N_hospitals] v;
    v = diag_pre_multiply(sigmas ,L) * z_matrix;

   // Compose the Random Ward Intercept
   vector[N_wards] alpha_wards;
   for(i in 1:N_wards){
    alpha_wards[i] = ward_a_bar + ward_sigma * zj[i];
   }
}
// Model 
model{
    // Correlated hospital intercept and Dose effect Priors 
    a_bar  ~ normal(-0.4, 0.1);   // Hospital intercept mean 
    b_bar  ~ normal(0.8,  0.1);   // Dose Effect mean  
    sigmas ~ exponential(1);
    L      ~ lkj_corr_cholesky(2);
    to_vector(z_matrix) ~ std_normal();
  
    // Ward intercept Priors 
    ward_a_bar ~ normal(0, 1);
    ward_sigma ~ exponential(1);
    zj         ~ std_normal();

    // Fixed Effects Priors
    beta_comorbidities ~ normal(-0.35,  0.01);
    beta_age           ~ normal(-0.015, 0.01);

    // Model Likelihood 
    vector[N] pi;
    for(i in 1:N){
        pi[i] = inv_logit(
            (a_bar + v[1, hospital_id[i]]) // Hospital Intercept
            +
            alpha_wards[ward_id[i]]        // Ward Intercepts
            +
            (b_bar + v[2, hospital_id[i]]) * dose[i] // Dose Random Effect 
            + 
            beta_comorbidities * comorbidities[i]
            +
            beta_age           * age[i]
        );
    }
    // Prior Switch 
    if(prior_only == 0){
        // Obervation Model Sample from Bernoulli Distribution 
        recovery ~ bernoulli(pi);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] pi;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        pi[i] = inv_logit(
                (a_bar + v[1, hospital_id[i]]) // Hospital Intercept
                +
                alpha_wards[ward_id[i]]        // Ward Intercepts
                +
                (b_bar + v[2, hospital_id[i]]) * dose[i] // Dose Random Effect 
                + 
                beta_comorbidities * comorbidities[i]
                +
                beta_age           * age[i]
        );
        log_lik[i] = bernoulli_lpmf(recovery[i] | pi[i]);
        y_rep[i]   = bernoulli_rng(pi[i]);
    }
    // Correlation Recovery Specific for the Correlated Models 
    matrix[2,2] Rho;
    Rho = multiply_lower_tri_self_transpose(L); 
}