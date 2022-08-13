// Jan 15, 2022

//July 4 2022 modification - adding ncp for both the slope and the intercept
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
  //int species[N];
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
  //real mu_grand;
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
  vector[Nspp] b_tilde; 
  real b_z;
  vector[Nspp] a_tilde; // intercept
  real a_z; // grand mean
  // vector[Nspp] a_;
  // vector[Nspp] b_;
	}
	
// transformed parameters {
//   vector[Nspp] a = a_z + sqrt(lam_interceptsa) * sigma_interceptsa * (L * a_); // implies: a ~ multi_normal
//   vector[Nspp] b = b_z + sqrt(lam_interceptsb) * sigma_interceptsb * (L * b_); // implies: b ~ multi_normal
// }
transformed parameters{
  vector[Nspp] a;
  vector[Nspp] b;
  // These braces allow us to define vcv_a locally to the `transformed parameters` block so that it doesn’t get saved in the Stan output
    //  Ideally vcv_a would be called something more descriptive like L_vcv_a
{
 matrix[Nspp,Nspp] L_vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
a = L_vcv_a * a_tilde;
  }
  
  {
 matrix[Nspp,Nspp] L_vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
b = L_vcv_b * b_tilde;
  }
}


model {
  real yhat[N];
  //matrix[Nspp,Nspp] L_vcv_a;     // phylogeny
  //matrix[Nspp,Nspp] L_vcv_b;     // phylogeny
  
  for(i in 1:N){
    yhat[i] = 
      a_z + a[species[i]] + b_z + b[species[i]] * year[i]; //astudy[study[i]] +
			     	}
  //L_vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  //L_vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
  a_tilde ~ normal(0,1);
  b_tilde ~ normal(0,1);
  
   ypred ~ normal(yhat, sigma_y);
  //astudy ~normal(0, astudy_sigma);

 // Prior
  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}









