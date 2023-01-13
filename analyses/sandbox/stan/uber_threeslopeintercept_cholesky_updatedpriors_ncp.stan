// Started January 12 by Deirdre
// the purpose of this code is to test the new phylogeny model on the ospree phylogeny model:
// 1. comment out the creation of sigma_mat from the function (line 16)
// 2. creating the matrix for the intercept and each slope (84-87)
// 3. add cholesky decomposition of the matricies line 95-98
// 4. add use multi_normal_cholesky distribution and rep_vector to return a matrix of size n_sp consisting of copies of the vcv matrix (101-104)
// Updated 23 Sep 2022 by Lizzie with better priors //

functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
   // matrix[rows(vcv),cols(vcv)] sigma_mat;  
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
      return(quad_form_diag(local_vcv, rep_vector(sigma, rows(vcv))));
    //sigma_mat = diag_matrix(rep_vector(sigma, rows(vcv)));
    //return(sigma_mat * local_vcv * sigma_mat);
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] y; 		// response
  vector[N] x1; 	// predictor (springtemp)
  matrix[n_sp,n_sp]Vphy;     // phylogeny
}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsbf;       
  real<lower=0> sigma_interceptsbf;   
  
  real b_zf; // 
  vector[n_sp] a; // intercept
  real a_z;
  
  real mu_b_force;
  vector[n_sp] b_force_ncp; // slope of forcing effect
  real <lower=0> sigma_b_force;

}

transformed parameters {
  vector[n_sp] b_force;
  
  b_force = mu_b_force + sigma_b_force*b_force_ncp;
}

model {
       real yhat[N];
       
        matrix[n_sp,n_sp] vcv_a;     // phylogeny
        matrix[n_sp,n_sp] vcv_bf;     // phylogeny

       
       	for(i in 1:N){
            yhat[i] = 
	      a[sp[i]] + b_force[sp[i]] * x1[i];
			     	}
			     	
	vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_bf = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsbf, sigma_interceptsbf));

  a ~ multi_normal_cholesky(rep_vector(a_z,n_sp), vcv_a); 
  b_force ~ multi_normal_cholesky(rep_vector(b_zf,n_sp), vcv_bf); 
  
  y ~ normal(yhat, sigma_y);

 // Priors -- keep in Stan code, better for reproducibility and runs faster
    a_z ~ normal(100, 10); // Same as before, seems okay
    b_zf ~ normal(-2, 10); // updated prior ... I think we should also try 0, 10


    // All below: same as before, seems okay
    lam_interceptsa ~ beta(1, 1);
    lam_interceptsbf ~ beta(1, 1);
  

    // I don't have a good sense of how to set these, so keeping a little wide
    sigma_interceptsa ~ normal(30, 20);
    sigma_interceptsbf ~ normal(1, 5);
    
    sigma_y ~ normal(10, 10); // updated prior
    
    //target += normal_lpdf(to_vector(b_force_ncp) | 0,1);
    b_force_ncp ~ normal(0,1);
    sigma_b_force ~ normal(0,10);
    
  
}


