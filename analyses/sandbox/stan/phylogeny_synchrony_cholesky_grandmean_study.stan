
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
  int<lower=1> Nspp;

  int<lower=1, upper= Nspp > species[N];

  int<lower=0> Nstudy; 
  int study[N]; 
  
  vector[N] ypred; 		// response
  vector[N] year; 	// predictor (year)
  matrix[Nspp, Nspp] Vphy;     // phylogeny
  // Priors
  
  real lam_interceptsa_prior_alpha;
  real lam_interceptsa_prior_beta;
  real sigma_interceptsa_prior_mu;
  real sigma_interceptsa_prior_sigma;
  real b_z_prior_mu;
  real b_z_prior_sigma;
  real lam_interceptsb_prior_alpha;
  real lam_interceptsb_prior_beta;
  real sigma_interceptsb_prior_mu;
  real sigma_interceptsb_prior_sigma;
  real sigma_y_prior_mu;
  real sigma_y_prior_sigma;  	
}

parameters {
  real mu_grand;
  
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
  vector[Nspp] b; // slope
  real b_z;
  vector[Nspp] a; // intercept

  real<lower=0> sigma_study; // variation in int among sp  
  real mu_study[Nstudy];
	}

model {
  real yhat[N];
  matrix[Nspp,Nspp] vcv_a;     // phylogeny
  matrix[Nspp,Nspp] vcv_b;     // phylogeny
  
  for(i in 1:N){
    yhat[i] = 
      mu_grand + a[species[i]] + mu_study[study[i]] + b[species[i]] * year[i]; //astudy[study[i]] +
			     	}
  vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
  a ~ multi_normal_cholesky(rep_vector(0,Nspp), vcv_a);
  b ~ multi_normal_cholesky(rep_vector(b_z, Nspp), vcv_b); 
  
   ypred ~ normal(yhat, sigma_y);

 // Prior
  mu_grand ~ normal (150, 50); 
  
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);
  
  sigma_study ~ normal(0,150);
  mu_study ~ normal(0, sigma_study);

}









