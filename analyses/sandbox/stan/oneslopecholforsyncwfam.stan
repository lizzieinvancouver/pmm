// Started Jan 15, 2022 by D. Loughnan 
// Updated 15 May 2022 by Lizzie 

/* Model identical to what Lizzie with help from W. Pearse wrote in 2020 (see ubermini_2.stan) but with ...
1) Priors added in from R 
2) Added lines to cholesky_decompose vcv (a and b)
3) Changed from multi_normal to multi_normal_cholesky
4) Geoff Legault removed sigma_mat from f(x), and just used sigma 
Some/all of above 2-4 changes seemed to have sped things up */

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
  int<lower=1, upper= Nfam > family[N];
  vector[N] ypred; 		// response
  vector[N] year; 	// predictor (year)
  matrix[Nfam, Nfam] Vphy;     // phylogeny
  
  // Priors
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
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
  vector[Nfam] b; // slope of forcing effect
  real b_z;
	}


model {
  real yhat[N];
  matrix[Nfam, Nfam] vcv_b;     // phylogeny
  
  for(i in 1:N){
    yhat[i] = 
      b[family[i]] * year[i]; 
			     	}
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
  b ~ multi_normal_cholesky(rep_vector(b_z, Nfam), vcv_b); 
  
   ypred ~ normal(yhat, sigma_y);

 // Prior
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}









