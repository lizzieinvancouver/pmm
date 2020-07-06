// uber pgls
// pretty sure this is wrong.... and it's not working

data {
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] x; 	// predictor
        matrix[n_sp,n_sp]Vphy;     // phylogeny
		
	}

parameters {
  matrix[N, N] sigma_y;    
  real<lower=0,upper=100> null_interceptsb;       
  real<lower=0,upper=100> lam_interceptsb;    
  real a; // intercept
  real b_force; // slope of forcing effect 
  real yhat[N];
	}

model {
       sigma_y = diag_matrix(rep_vector(null_interceptsb, n_sp)) + lam_interceptsb*Vphy; 
       yhat = a + b_force * x;
        
	a ~ normal(0,10);
        b_force ~ normal(0,3);
        null_interceptsb ~ normal(0, 20);
	lam_interceptsb ~ normal(0, 20);
	y ~ multi_normal(yhat, sigma_y);

}

