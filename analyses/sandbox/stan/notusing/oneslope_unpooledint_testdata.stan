// started January 22, 2022 

//purpose of this code is to model phylogeny on the slopes and unpooled intercepts
//code adapted from uber_onleslope.stan - but with unpooled spp intercepts 

functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
    matrix[rows(vcv),cols(vcv)] sigma_mat;  
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
    sigma_mat = diag_matrix(rep_vector(sigma, rows(vcv)));
    return(sigma_mat * local_vcv * sigma_mat);
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] y; 		// response
  vector[N] x1; 	// predictor
  matrix[n_sp,n_sp]Vphy;     // phylogeny
  // Priors
  // real mu_a_prior_mu;
  // real mu_a_prior_sigma;
  real a_prior_sigma;
  // real sigma_a_prior_sigma;
  real b_z_prior_mu;
  real b_z_prior_sigma;
  real lam_interceptsb_prior_alpha;
  real lam_interceptsb_prior_beta;
  real sigma_interceptsb_prior_mu;
  real sigma_interceptsb_prior_sigma;
  real sigma_y_mu_prior;
  real sigma_y_mu_sigma;  
	}

parameters {
  
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;    
  vector[n_sp] b_force; // slope of forcing effect
  real b_z;
  
  vector[n_sp] a;
  //real mu_a;
  //real <lower = 0> sigma_a;
  	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = a[sp[i]] + b_force[sp[i]] * x1[i];
	}
  
  b_force ~ multi_normal(rep_vector(b_z,n_sp), lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); 
  y ~ normal(yhat, sigma_y);
  
  // Priors
  a ~ normal(0, a_prior_sigma); 
  // mu_a ~ normal(mu_a_prior_mu, mu_a_prior_sigma);
  // sigma_a ~ normal(sigma_a_prior_mu, sigma_a_prior_sigma);
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sigma_y ~ normal(sigma_y_mu_prior, sigma_y_mu_sigma);
  
}

