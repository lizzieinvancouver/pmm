// Started Jan 15, 2022 by D. Loughnan 
// Updated 15 May 2022 by Lizzie 

/* Model identical to what Lizzie with help from W. Pearse wrote in 2020 (see ubermini_2.stan) but with ...
1) Priors added in from R 
2) Added lines to cholesky_decompose vcv (a and b)
3) Changed from multi_normal to multi_normal_cholesky
4) Geoff Legault removed sigma_mat from f(x), and just used sigma 
Some/all of above 2-4 changes seemed to have sped things up */

/* Notes and sample code here from Mike Betancourt */
/* This code DOES NOT run, it's not complete. */

functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){ //input: vcv matrix (correl), lam, sig which are sampled from dist
    matrix[rows(vcv),cols(vcv)] local_vcv; //blank matrix - size vcv
    local_vcv = vcv * lambda; // taking the correl matrix and multi by lambda - scales correl matrix; how far or close evo are things rel to phylo dist
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i]; //ensure dist on diag are rescaled to 1; undoing the multi of diagonals to set back to 1 (i,i)
   return(quad_form_diag(local_vcv, rep_vector(sigma, rows(vcv)))); //if know correl mat and sigma - could build covariance matrix, really only build two sep matrices to calculate this
  }
}

data {
  int<lower=1> N;
  int<lower=1> Nspp;
  int<lower=1, upper= Nspp > species[N];
  vector[N] ypred; 		// response
  vector[N] year; 	// predictor (year)
  matrix[Nspp, Nspp] Vphy;     // phylogeny
  
  // Priors
  real a_z_prior_mu; 
  real a_z_prior_sigma;
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
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
  vector[Nspp] b; // slope of forcing effect
  real b_z;
  vector[Nspp] a; // intercept
  real a_z;
	}


model {
  real yhat[N];
  matrix[Nspp,Nspp] vcv_a;     // phylogeny
  matrix[Nspp,Nspp] vcv_b;     // phylogeny
  /* START of Mike adding notes on how to make the phylogeny model have another level
  ... but needs edits above (made in other files: oneinterceptcholforsyncfam.stan for now) to run. */
  for(i in 1:N){
    yhat[i] = 
      a_z + a[family[i]] + a_species[species[i]] + b[family[i]] * year[i]; 
			     	}
  vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
  afamily ~ multi_normal_cholesky(rep_vector(0, Nspp), vcv_a);
  bfamily ~ multi_normal_cholesky(rep_vector(0, Nspp), vcv_b); 
  
  for (s in 1:N_species) {
    a_species[s] ~ normal(0, tau_family[family[s]])
  }
  
  // Demo Start (by Mike of how to vectorize code)
  int<lower=1> N;  // Number of observations
  int<lower=1> K;  // Number of contexts
  int<lower=1, upper=K> context[N]; // Observation to context mapping
  
  real y[N];
  real theta[K];
  
  for (n in 1:N) {
    y[n] ~ normal(theta[context[n]], sigma);
  }
  
  theta[context[n]] -> real
  theta[context]    -> array[N] real
  
  y ~ normal(theta[context], sigma)
  
  y ~ normal(alpha[context] + beta[context] .* year)
  
  
  y ~ normal(alpha_family[family] + alpha_species[species], sigma);
  y ~ normal(alpha_family[family[species]] + alpha_species[species], sigma);
  
  // Demo End
  
   ypred ~ normal(yhat, sigma_y);

 // Prior

  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}









