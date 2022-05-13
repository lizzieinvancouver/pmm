// Jan 15, 2022
// Given that the NCP model that Ben Goodrich suggested ran faster on real data than with the fake data, I am testing whether Geoff's new code, which takes ~10 days to run on the test data also runs faster
//this model is not ncp, Geoff also modified it by removing the sigma_mat from the function

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
  //int<lower=1> Nstudy;
  //int species[N];
  //int study[N];
  int<lower=1, upper= Nspp > species[N];
  //int<lower=1, upper= Nstudy > study[N];
  vector[N] ypred; 		// response
  vector[N] year; 	// predictor (year)
  matrix[Nspp, Nspp] Vphy;     // phylogeny
  // Priors
  
  real a_z_prior_mu; // mu_grand species intercept - w/o study
  real a_z_prior_sigma;
  //real astudy_sigma_prior_mu;
  //real astudy_sigma_prior_sigma;
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
  //real a_prior_mu;
  //real a_prior_sigma;
  //real b_prior_mu;
  //real b_prior_sigma;
}

// transformed data {
//   matrix[Nspp, Nspp] L = cholesky_decompose(Vphy);
// }

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
  // vector[Nspp] a_;
  // vector[Nspp] b_;
  //vector[Nstudy] astudy; // intercept
  //real astudy_sigma; //study level effect offset from multinormal; not grand mean
	}
	
// transformed parameters {
//   vector[Nspp] a = a_z + sqrt(lam_interceptsa) * sigma_interceptsa * (L * a_); // implies: a ~ multi_normal
//   vector[Nspp] b = b_z + sqrt(lam_interceptsb) * sigma_interceptsb * (L * b_); // implies: b ~ multi_normal
// }

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
  //astudy ~normal(0, astudy_sigma);

 // Prior

  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  //astudy_sigma ~ normal(astudy_sigma_prior_mu, astudy_sigma_prior_sigma);
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}









