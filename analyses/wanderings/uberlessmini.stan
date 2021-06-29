// uber less mini test
// Phylogeny on the intercept

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
  real<lower=0,upper=100> null_interceptsb;       
  real<lower=0,upper=100> lam_interceptsb;    
  vector[n_sp] a; // intercepts
  real b_force; // slope of forcing effect 
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
	a[sp[i]] + b_force * x[i];
			     	}

	a ~ multi_normal(rep_vector(0,n_sp), diag_matrix(rep_vector(null_interceptsb, n_sp)) + lam_interceptsb*Vphy); 

        null_interceptsb ~ normal(0, 20);
	lam_interceptsb ~ normal(0, 20);
        b_force ~ normal(0,3);
	y ~ normal(yhat, sigma_y);

}

