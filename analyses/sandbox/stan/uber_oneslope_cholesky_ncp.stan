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
  vector[N] x1; 	// predictor (year)
  
  int<lower=1> Nspp;
  int<lower=1, upper= Nspp > species[N];
  matrix[Nspp, Nspp] Vphy;     // phylogeny
}


parameters {
  //real mu_grand;
  real<lower=0> sigma_y;    
  
  real a_z; // grand mean
  vector[Nspp] a; // intercept

  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  
  real b_z;
  vector[Nspp] b_tilde; 
  
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;   
	}

transformed parameters {
  // Add back in sigma scaling
                 
  vector[Nspp] b =   sigma_interceptsb
                   * cholesky_decompose(unscaled_lambda_vcv(Vphy, lam_interceptsb)) * b_tilde;        
}
	
model {
  vector[N] yhat;
       
        matrix[Nspp,Nspp] vcv_a;     // phylogeny
     // phylogeny

       
       	for(i in 1:N){
            yhat[i] = 
	      a[species[i]] + b[species[i]] * x1[i];
			     	}

	      
//  vector[N] yhat = a_z + a[species] + ((b_z + b[species])).* x1;
  vcv_a = cholesky_decompose(unscaled_lambda_vcv(Vphy, lam_interceptsa));
  a ~ multi_normal_cholesky(rep_vector(a_z,Nspp), vcv_a); 
  
  a_z ~ normal(100, 10); // Same as before, seems okay
  
  lam_interceptsa ~ beta(1,1);
  sigma_interceptsa ~ normal(30,20);
  
  b_z ~ normal(-2,10);
  b_tilde ~ normal(0,1);
  
  lam_interceptsb ~ beta(1,1);
  sigma_interceptsb ~ normal(1,5);
  
  sigma_y ~ normal(10,10);

  yobs ~ normal(yhat, sigma_y);
  
 
}









