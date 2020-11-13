// uber mini test

functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
    return(local_vcv * sigma);
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
  
  matrix[1, n_sp] z;       //matrix oz z values used in ncp Cholesky decomposition 
	}
	
transformed parameters{
  //decomposing the matrix made using fuinction lambda_vcv
  cholesky_factor_corr[n_sp] L_lambda_vcv; 
  L_lambda_vcv = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));

}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
		b_force[sp[i]] * x[i];
			     	}
  // this distribution used teh decomposed matrix (ncp happens, but we dont see the z value with this function)
	b_force ~ multi_normal_cholesky(rep_vector(b_z,n_sp), L_lambda_vcv); 

        sigma_interceptsb ~ normal(0, 1);
	lam_interceptsb ~ normal(0, 1);
	b_z ~ normal(0, 0.5);
	y ~ normal(yhat, sigma_y);

}

