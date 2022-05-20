// Started 18 May 2022 by Lizzie 

/* Based on oneslopeinterceptcholforsync.stan but ...
1) made it intercept-only and pulled out grand mean 
2) Tried to add the family - species hierachical levels */

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
  int<lower=1> Nfam;
  int<lower=1> Nspp;
  int<lower=1, upper= Nspp > species[N];
  int<lower=1, upper= Nfam > family[N];
  vector[N] ypred; 		// response
  matrix[Nfam, Nfam] Vphy;     // phylogeny
  
  // Priors
  real a_z_prior_mu; 
  real a_z_prior_sigma;
  real lam_interceptsa_prior_alpha;
  real lam_interceptsa_prior_beta;
  real sigma_interceptsa_prior_mu;
  real sigma_interceptsa_prior_sigma;
  real sigma_y_prior_mu;
  real sigma_y_prior_sigma;  	
}



parameters {
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  vector[Nspp] a_species_raw; // 
  vector[Nfam] a; // intercept
  vector[Nfam] tau_family; 
  real a_z;
	}

transformed parameters {
  vector[Nspp] a_species; // 
  for (s in 1:Nspp) {
     a_species[s] = tau_family[family[s]] * a_species_raw[s]; /* NCP attempt*/
  }
}

model {
  real yhat[N];
  matrix[Nfam,Nfam] vcv_a;     // phylogeny

  for(i in 1:N){
    yhat[i] = 
      a_z + a[family[i]] + a_species[species[i]]; // Yes, this should include both family and species
			     	}
   vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
   a ~ multi_normal_cholesky(rep_vector(0, Nfam), vcv_a); // family level 
   a_species_raw ~ normal(0, 1); 
    
   ypred ~ normal(yhat, sigma_y);

 // Prior

  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}
