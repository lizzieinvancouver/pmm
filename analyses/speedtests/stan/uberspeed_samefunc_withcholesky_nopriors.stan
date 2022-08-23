/* Copied ubermini_2_biggerpriors.stan to start, then edited
By Lizzie and Dinara Mamatova
uber mini test with phylogenetic partial pooling on slope *and* intercept, so:
  ypred = a[phylo] + b_force[phylo]*x
y ~ normal(ypred, sigma_y) */
  
  functions {
    matrix lambda_vcv(matrix vcv, real lambda, real sigma){
      matrix[rows(vcv),cols(vcv)] local_vcv;
      matrix[rows(vcv),cols(vcv)] sigma_mat;
      local_vcv = vcv * lambda;
      for(i in 1:rows(local_vcv))
        local_vcv[i,i] = vcv[i,i];
      sigma_mat = diag_matrix(rep_vector(sigma, rows(vcv)));
      return(sigma_mat * local_vcv * sigma_mat);
      // Another way of writing the output of function without saving sigma_mat variable
      // return(quad_form_diag(local_vcv, rep_vector(sigma, rows(vcv))));
    }
  }

data {
  int<lower=1> N;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] y; 		// response
  vector[N] x; 	// predictor
  matrix[n_sp,n_sp]Vphy;     // phylogeny
  
  // Priors
  // real a_z_prior_mu;
  // real a_z_prior_sigma;
  // real lam_interceptsa_prior_alpha;
  // real lam_interceptsa_prior_beta;
  // real sigma_interceptsa_prior_mu;
  // real sigma_interceptsa_prior_sigma;
  // real b_z_prior_mu;
  // real b_z_prior_sigma;
  // real lam_interceptsb_prior_alpha;
  // real lam_interceptsb_prior_beta;
  // real sigma_interceptsb_prior_mu;
  // real sigma_interceptsb_prior_sigma;
  // real sigma_y_mu_prior;
  // real sigma_y_mu_sigma;
}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;    
  vector[n_sp] b_force; // slope of forcing effect
  real b_z;
  vector[n_sp] a; // intercept
  real a_z;
  real<lower=0> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;    
}

model {
  
  real yhat[N];
  
  // Phylogeny for Cholesky decomposition
  matrix[n_sp,n_sp] vcv_a;
  matrix[n_sp,n_sp] vcv_b;
  
  for(i in 1:N){
    yhat[i] = a[sp[i]] +
      b_force[sp[i]] * x[i];
  }
  /* print(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); */
    
    // Use Cholesky decomposition
  vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));

  a ~ multi_normal_cholesky(rep_vector(a_z,n_sp), vcv_a);
  b_force ~ multi_normal_cholesky(rep_vector(b_z, n_sp), vcv_b);
  
  // b_force ~ multi_normal(rep_vector(b_z,n_sp), lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
  // a ~ multi_normal(rep_vector(a_z,n_sp), lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  
  y ~ normal(yhat, sigma_y);
  
  // Generate paramaters using constant values
  sigma_interceptsb ~ normal(0, 1);
  lam_interceptsb ~ beta(2, 2);
  b_z ~ normal(0, 5);
  sigma_interceptsa ~ normal(0, 1);
  lam_interceptsa ~ beta(2, 2);
  a_z ~ normal(0, 5);
  sigma_y ~ normal(0, 0.5);
  
  // Generate parameters using priors
  // sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  // lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  // b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  // sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  // lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  // a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  // sigma_y ~ normal(sigma_y_mu_prior, sigma_y_mu_sigma);
  
}

