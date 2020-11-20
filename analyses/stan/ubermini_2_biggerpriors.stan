// uber mini test

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
	vector[N] x; 	// predictor
        matrix[n_sp,n_sp]Vphy;     // phylogeny
		
	}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;    
  vector[n_sp] b_force; // slope of forcing effect
  real b_z;
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
		b_force[sp[i]] * x[i];
			     	}
  /* print(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); */
  
  b_force ~ multi_normal(rep_vector(b_z,n_sp), lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); 

  sigma_interceptsb ~ normal(0, 1);
  lam_interceptsb ~ beta(2, 2);
  b_z ~ normal(0, 5);
  y ~ normal(yhat, sigma_y);
  sigma_y ~ normal(0, .5);

}

