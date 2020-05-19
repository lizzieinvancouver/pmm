// OSPREE analysis
// Simplified version of nointe_2level.stan 
// with phylogeny on intercept and turned off partial pooling on the slopes

// NEED TO keep working on this ...
// See https://groups.google.com/forum/#!topic/stan-users/Irv9RWDCpQE

data {
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] force; 	// predictor
        matrix[n_sp,n_sp]Vphy;     // phylogeny
		
	}

parameters {
  // real mu_a_sp;   
  //real mu_b_force_sp;   
  //real<lower=0> sigma_a_sp; 
  // real<lower=0> sigma_b_force_sp; 
  real<lower=0> sigma_y; 
  real<lower=0,upper=100> null_intercepts;       
  real<lower=0,upper=100> lam_intercepts;       
  vector[n_sp] a_sp; // intercept for species
  real b_force; // slope of forcing effect 
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = a_sp[sp[i]] + // indexed with species
		b_force * force[i];
			     	}

	a_sp ~ multi_normal(rep_vector(0,n_sp), diag_matrix(rep_vector(null_intercepts, n_sp)) + lam_intercepts*Vphy); 

        //mu_b_force_sp ~ normal(0, 50);
        // sigma_b_force_sp ~ normal(0, 10);
   
        //mu_a_sp ~ normal(0, 50);
        // sigma_a_sp ~ normal(0, 10);
	b_force ~ normal(0, 10);
	
	y ~ normal(yhat, sigma_y);

}


