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
	vector[N] x1; 	// predictor (forcing)
	vector[N] x2; 	// predictor (chilling)
	vector[N] x3; 	// predictor (photo)
        matrix[n_sp,n_sp]Vphy;     // phylogeny
		
	}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0> lam_interceptsbf;       
  real<lower=0> sigma_interceptsbf;   
  real<lower=0> lam_interceptsbc;       
  real<lower=0> sigma_interceptsbc; 
  real<lower=0> lam_interceptsbp;       
  real<lower=0> sigma_interceptsbp;  
  vector[n_sp] b_force; // slope of forcing effect
  real b_zf;
  vector[n_sp] b_chill; // slope of chilling effect
  real b_zc;
  vector[n_sp] b_photo; // slope of photo effect
  real b_zp;
  vector[n_sp] a; // intercept
  real a_z;
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
		a[sp[i]] + b_force[sp[i]] * x1[i] + 
                b_chill[sp[i]] * x2[i] + b_photo[sp[i]] * x3[i]; //
			     	}
  /* print(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb)); */
  a ~ multi_normal(rep_vector(a_z,n_sp), lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa)); 
  b_force ~ multi_normal(rep_vector(b_zf,n_sp), lambda_vcv(Vphy, lam_interceptsbf, sigma_interceptsbf)); 
  b_chill ~ multi_normal(rep_vector(b_zc,n_sp), lambda_vcv(Vphy, lam_interceptsbc, sigma_interceptsbc)); 
  b_photo ~ multi_normal(rep_vector(b_zp,n_sp), lambda_vcv(Vphy, lam_interceptsbp, sigma_interceptsbp)); 
  y ~ normal(yhat, sigma_y);

 // Priors
  a_z ~ normal(4, 1);
  lam_interceptsa ~ beta(2, 2);
  sigma_interceptsa ~ normal(0.2, 0.1);
  b_zf ~ normal(0.6, 0.1);
  lam_interceptsbf ~ beta(2, 2);
  sigma_interceptsbf ~ normal(0.1, 0.05);
  b_zc ~ normal(0.88, 0.1);
  lam_interceptsbc ~ beta(5, 2);
  sigma_interceptsbc ~ normal(0.07, 0.02);
  b_zp ~ normal(1.1, 0.4);
  lam_interceptsbp ~ beta(2, 2);
  sigma_interceptsbp ~ normal(0.3, 0.05);
  sigma_y ~ normal(0.02, 0.01);

}
