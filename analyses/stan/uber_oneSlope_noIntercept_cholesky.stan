// January 25: turning off the phylogeny on the intercept to see if it speeds things up
functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
    return(quad_form_diag(local_vcv, rep_vector(sigma, rows(vcv))));
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] y; 		// response
  vector[N] x1; 	// predictor
  matrix[n_sp,n_sp] Vphy;     // phylogeny
  // Priors
  real mu_a_prior_mu;
  real mu_a_prior_sigma;
  real sigma_a_prior_mu;
  real sigma_a_prior_sigma;
  // real a_z_prior_mu;
  // real a_z_prior_sigma;
  // real lam_interceptsa_prior_alpha;
  // real lam_interceptsa_prior_beta;
  // real sigma_interceptsa_prior_mu;
  // real sigma_interceptsa_prior_sigma;
  real b_zf_prior_mu;
  real b_zf_prior_sigma;
  real lam_interceptsbf_prior_alpha;
  real lam_interceptsbf_prior_beta;
  real sigma_interceptsbf_prior_mu;
  real sigma_interceptsbf_prior_sigma;
  real sigma_y_prior_mu;
  real sigma_y_prior_sigma;  

} 

parameters {
  real<lower=0> sigma_y;    
  // real<lower=0, upper=1> lam_interceptsa;       
  // real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsbf;       
  real<lower=0> sigma_interceptsbf;    
  vector[n_sp] b_force; // slope of forcing effect
  real b_zf;
  vector[n_sp] a; // intercept
  real mu_a;
  real<lower = 0> sigma_a; 
	}

model {
  real yhat[N];
  // matrix[n_sp,n_sp] vcv_a;     // phylogeny
  matrix[n_sp,n_sp] vcv_b;     // phylogeny

  for(i in 1:N){
    yhat[i] = 
      a[sp[i]] + b_force[sp[i]] * x1[i];
  }

  // vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsbf, sigma_interceptsbf));
 
  // a ~ multi_normal_cholesky(rep_vector(a_z,n_sp), vcv_a);
  a ~ normal(mu_a, sigma_a);
  b_force ~ multi_normal_cholesky(rep_vector(b_zf, n_sp), vcv_b);

  y ~ normal(yhat, sigma_y);
  
  // Priors
  mu_a ~ normal(mu_a_prior_mu, mu_a_prior_sigma);
  sigma_a ~ normal(sigma_a_prior_mu, sigma_a_prior_sigma);
  // lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  // sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  b_zf ~ normal(b_zf_prior_mu, b_zf_prior_sigma);
  lam_interceptsbf ~ beta(lam_interceptsbf_prior_alpha, lam_interceptsbf_prior_beta);
  sigma_interceptsbf ~ normal(sigma_interceptsbf_prior_mu, sigma_interceptsbf_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}

