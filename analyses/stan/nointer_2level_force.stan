// OSPREE analysis
// Simplified version of nointe_2level.stan

data {
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] force; 	// predictor
		
	}

parameters {
  real mu_a_sp;   
  real mu_b_force_sp;   
  real<lower=0> sigma_a_sp; 
  real<lower=0> sigma_b_force_sp; 
  real<lower=0> sigma_y; 

  real a_sp[n_sp]; // intercept for species
  real b_force[n_sp]; // slope of forcing effect 
	}

transformed parameters {
   real yhat[N];
       	for(i in 1:N){
            yhat[i] = a_sp[sp[i]] + // indexed with species
		b_force[sp[i]] * force[i];
			     	}

	}

model {

	a_sp ~ normal(mu_a_sp, sigma_a_sp); 
	b_force ~ normal(mu_b_force_sp, sigma_b_force_sp); 

        mu_b_force_sp ~ normal(0, 50);
        sigma_b_force_sp ~ normal(0, 10);
   
        mu_a_sp ~ normal(0, 50);
        sigma_a_sp ~ normal(0, 10);
	//b_force ~ normal(0, 10);
	
	y ~ normal(yhat, sigma_y);

}

generated quantities{
   real y_ppc[N];
   for (n in 1:N)
      y_ppc[n] = a_sp[sp[n]] + 
		b_force[sp[n]] * force[n];
    for (n in 1:N)
      y_ppc[n] = normal_rng(y_ppc[n], sigma_y);

}
