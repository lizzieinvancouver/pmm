// Jan 15, 2022

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
  
  real a_z_prior_mu; // mu_grand species intercept - w/o study
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
  vector[Nspp] b; // slope 
  real b_z;
  vector[Nspp] a; // intercept
  real a_z; 
	}

model {
  real yhat[N];
  matrix[Nspp,Nspp] vcv_a;     // phylogeny
  matrix[Nspp,Nspp] vcv_b;     // phylogeny
  
  for(i in 1:N){
    yhat[i] = 
      a[species[i]] + b[species[i]] * year[i]; //astudy[study[i]] +
			     	}
  vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
  a ~ multi_normal_cholesky(rep_vector(a_z,Nspp), vcv_a);
  b ~ multi_normal_cholesky(rep_vector(b_z, Nspp), vcv_b); 
  
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









