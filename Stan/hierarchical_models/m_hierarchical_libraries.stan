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
// Transform data 
transformed data {
   // Scale numerical covariates 
   vector[N] z_rain_mm;
   vector[N] z_held_book;
   vector[N] z_staff_count;
   z_rain_mm = (rain_mm - mean(rain_mm)) / sd(rain_mm);
   z_held_book = (held_book - mean(held_book)) / sd(held_book);
   z_staff_count = (staff_count - mean(staff_count)) / sd(staff_count);
}
// Model Paramters 
parameters{
    // Random Intercept 
    real z_mean_book_loans;
    real <lower = 0.001> sigma_book_count;
    vector[n_libraries] zj;

    // Fixed covariates effects 
    real c_beta_holiday;
    real z_beta_rain;
    real z_beta_held;
    real z_beta_staff;
}
// Transform NCP of the model Random Intercept 
transformed parameters {
   vector[n_libraries] alpha_j;
   
   for(i in 1:n_libraries){
    alpha_j[i] = z_mean_book_loans + sigma_book_count * zj[i];
   }
}
// Model 
model{
    // Random Intercept Priors 
    z_mean_book_loans  ~ normal(log(140), 2);
    sigma_book_count ~ exponential(1);
    zj               ~ std_normal();

    // Fixed covariates effects Priors 
    c_beta_holiday ~ normal(0, 1);
    z_beta_rain    ~ normal(0, 1);
    z_beta_held    ~ normal(0, 1);
    z_beta_staff   ~ normal(0, 1);

    // Model Likelihood 
    if(prior_only == 0){
        vector[N] lambda;
        for(i in 1:N){
            lambda[i] = exp(
                alpha_j[library_id[i]] 
                + c_beta_holiday * is_school_holiday[i]
                + z_beta_rain    * z_rain_mm[i]
                + z_beta_held    * z_held_book[i]
                + z_beta_staff   * z_staff_count[i]
            );
        }
    // Observation Poisson model 
    weekly_book_loans ~ poisson(lambda);
    }
}
// Minimal Generated Quantites 
generated quantities {
    // natural-scale quantities
    real c_beta_rain  = z_beta_rain  / sd(rain_mm);
    real c_beta_held  = z_beta_held / sd(held_book);
    real c_beta_staff = z_beta_staff / sd(staff_count);
    real mean_book_loans = z_mean_book_loans - c_beta_rain  * mean(rain_mm) - c_beta_held  * mean(held_book) - c_beta_staff * mean(staff_count);

    // Replications 
    vector[N] lambda;
    vector[N] log_lik;
    vector[N] y_rep;
    for(i in 1:N){
        lambda[i] = exp(
                alpha_j[library_id[i]] 
                + c_beta_holiday * is_school_holiday[i]
                + z_beta_rain    * z_rain_mm[i]
                + z_beta_held    * z_held_book[i]
                + z_beta_staff   * z_staff_count[i]
        );
        log_lik[i] = poisson_lpmf(weekly_book_loans[i] | lambda[i]);
        y_rep[i]   = poisson_rng(lambda[i]);
    }
}