// uber mini test

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
  vector[n_sp] b_force; // slope of forcing effect 
	}

model {
       real yhat[N];
       	for(i in 1:N){
            yhat[i] = 
		b_force[sp[i]] * x[i];
			     	}

	b_force ~ multi_normal(rep_vector(0,n_sp), diag_matrix(rep_vector(null_interceptsb, n_sp)) + lam_interceptsb*Vphy); 

        null_interceptsb ~ normal(0, 20);
	lam_interceptsb ~ normal(0, 20);
	y ~ normal(yhat, sigma_y);

}

