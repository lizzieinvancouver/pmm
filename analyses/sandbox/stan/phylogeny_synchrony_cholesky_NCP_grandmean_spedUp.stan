// Jan 15, 2022

functions {
  // BETAN: This the old lambda_vcv but with sigma factored out
  // input: vcv matrix (Phylogenetic correlation matrix)
  // input: lam (???)
   matrix unscaled_lambda_vcv(matrix vcv, real lambda){
    // taking the correl matrix and multi by lambda - scales correl matrix;
    // how far or close evo are things rel to phylo dist
    matrix[rows(vcv), cols(vcv)] local_vcv = vcv * lambda;

    //ensure dist on diag are rescaled to 1; undoing the multi of diagonals
    // to set back to 1 (i,i)
    for (i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
      return(local_vcv);
  }
}

data {
  int<lower=1> N;
  vector[N] yobs; 		// response
  vector[N] year; 	// predictor (year)
  
  int<lower=1> Nspp;
  int<lower=1, upper= Nspp > species[N];
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
  
  real a_z; // grand mean
  vector[Nspp] a_tilde; // intercept

  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  
  real b_z;
  vector[Nspp] b_tilde; 
  
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
	}

transformed parameters {
  // Add back in sigma scaling
  vector[Nspp] a =   sigma_interceptsa
                   * cholesky_decompose(unscaled_lambda_vcv(Vphy, lam_interceptsa)) * a_tilde;
                   
  vector[Nspp] b =   sigma_interceptsb
                   * cholesky_decompose(unscaled_lambda_vcv(Vphy, lam_interceptsb)) * b_tilde;        
}
	
model {
  vector[N] mu = a_z + a[species] + ((b_z + b[species])).* year;

  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  a_tilde ~ normal(0, 1);
  
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  b_tilde ~ normal(0,1);
  
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  
  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

  yobs ~ normal(mu, sigma_y);
  
  //real yhat[N];
  //matrix[Nspp,Nspp] vcv_a;     // phylogeny
 //matrix[Nspp,Nspp] vcv_b;     // phylogeny
  
//   for(i in 1:N){
//     yhat[i] = 
//       a_z + a[species[i]] + b_z + b[species[i]] * year[i]; //astudy[study[i]] +
// 			     	}
  //vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  //vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
 
 //  a_tilde ~ normal(0,1);
 //  b_tilde ~ normal(0,1);
 //  
 //   yobs ~ normal(yhat, sigma_y);
 //  //astudy ~normal(0, astudy_sigma);
 // 
 // // Prior
 //  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
 //  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
 //  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
 //  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
 //  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
 //  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
 //  sigma_y ~ normal(sigma_y_prior_mu, sigma_y_prior_sigma);

}









