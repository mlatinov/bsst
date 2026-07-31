// Input data Layer 
data{
    // Number of observations and prior switch 
    int <lower = 1> N;
    int <lower = 0, upper = 1> prior_only;

    // Covariates 
    vector[N] is_school_holiday;
    vector[N] rain_mm;
    vector[N] held_book;
    vector[N] staff_count;

    // Indexes 
    int <lower = 2> n_libraries;
    array[N] int library_id;

    // Outcome 
    array[N] int weekly_book_loans;
}
// Model Paramters 
parameters{
    // Random Intercept 
    real mean_book_loans;
    real <lower = 0.001> sigma_book_count;
    vector[n_libraries] zj;

    // Fixed covariates effects 
    real c_beta_holiday;
    real c_beta_rain;
    real c_beta_held;
    real c_beta_staff;
}
// Transform NCP of the model Random Intercept 
transformed parameters {
   vector[n_libraries] alpha_j;
   
   for(i in 1:n_libraries){
    alpha_j[i] = mean_book_loans + sigma_book_count * zj[i];
   }
}
// Model 
model{
    // Random Intercept Priors 
    mean_book_loans  ~ normal(log(140), 2);
    sigma_book_count ~ exponential(1);
    zj               ~ std_normal();

    // Fixed covariates effects Priors 
    c_beta_holiday ~ normal(0.18, 0.1);
    c_beta_rain    ~ normal(0.005, 0.01);
    c_beta_held    ~ normal(0.04, 0.1);
    c_beta_staff   ~ normal(0.06, 0.1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] lambda;
        for(i in 1:N){
            lambda[i] = exp(
                alpha_j[library_id[i]] 
                + c_beta_holiday * is_school_holiday[i]
                + c_beta_rain    * rain_mm[i]
                + c_beta_held    * held_book[i]
                + c_beta_staff   * staff_count[i]
            );
        }
    // Observation Poisson model 
    weekly_book_loans ~ poisson(lambda);
    }
}
// Minimal Generated Quantites 
generated quantities {
    vector[N] lambda;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        lambda[i] = exp(
                alpha_j[library_id[i]] 
                + c_beta_holiday * is_school_holiday[i]
                + c_beta_rain    * rain_mm[i]
                + c_beta_held    * held_book[i]
                + c_beta_staff   * staff_count[i]
        );
        log_lik[i] = poisson_lpmf(weekly_book_loans[i] | lambda[i]);
        y_rep[i]   = poisson_rng(lambda[i]);
    }
}