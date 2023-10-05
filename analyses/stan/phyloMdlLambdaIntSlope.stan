// Oct 4, 2023
// Started by D Loughnan
//Aim of this model --- to check whether model with lambda out performs model without
// This is the lambda model--- lambda on the slope and intercept---one cue (x)
// Testing different values lambda (lambda slope = lambda intercept): 0, 0.2,0.8


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
  vector[N] x1; 	// predictor (forcing)
  matrix[n_sp,n_sp]Vphy;     // phylogeny
}

parameters {
  real<lower=0> sigma_y;    
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
  vector[n_sp] b; // slope of forcing effect
  real b_z;

  vector[n_sp] a; // intercept
  real a_z;

}

model {
       real yhat[N];
       
        matrix[n_sp,n_sp] vcv_a;     // phylogeny
        matrix[n_sp,n_sp] vcv_b;     // phylogeny
       
       	for(i in 1:N){
            yhat[i] = 
	      a[sp[i]] + b[sp[i]] * x1[i];
			     	}
			     	
	vcv_a = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa));
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
  
  a ~ multi_normal_cholesky(rep_vector(a_z,n_sp), vcv_a); 
  b ~ multi_normal_cholesky(rep_vector(b_z,n_sp), vcv_b); 

  
  y ~ normal(yhat, sigma_y);

 // Priors -- keep in Stan code, better for reproducibility and runs faster
    a_z ~ normal(30, 10); // Same as before, seems okay
    b_z ~ normal(-2, 10); // updated prior ... I think we should also try 0, 10

    // All below: same as before, seems okay
    lam_interceptsa ~ beta(1, 1);
    lam_interceptsb ~ beta(1, 1);


    // I don't have a good sense of how to set these, so keeping a little wide
    sigma_interceptsa ~ normal(30, 20);
    sigma_interceptsb ~ normal(1, 5);
    
    sigma_y ~ normal(10, 10); // updated prior
  
}









