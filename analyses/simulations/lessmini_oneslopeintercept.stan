// based off ubermini_2_biggerpriors.stan
// but add an intercept ... 

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
		
	}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;    
  vector[n_sp] b_force; // slope of forcing effect
  real b_z;
  vector[n_sp] a; // intercept
  real a_z;
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
		a[sp[i]] + b_force[sp[i]] * x1[i];
			     	}
  a ~ multi_normal(rep_vector(a_z,n_sp), lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa)); 
  b_force ~ multi_normal(rep_vector(b_z,n_sp), lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); 
  y ~ normal(yhat, sigma_y);
  
  // Priors
  a_z ~ normal(4, 1);
  lam_interceptsa ~ beta(2, 2);
  sigma_interceptsa ~ normal(0.2, 0.1);
  b_z ~ normal(0.6, 0.1);
  lam_interceptsb ~ beta(2, 2);
  sigma_interceptsb ~ normal(0.1, 0.05);
  sigma_y ~ normal(0.02, 0.01);

}

