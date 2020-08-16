// uber pgls
// no longer sure what is wrong here ... compare to pgls_lemoine.stan

data {
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] x; 	// predictor
        matrix[n_sp, n_sp] Lmat; // (nsp by nsp and all 1s with 0s on diaganols) needed to make sure lambda is only calculated on off-diagonals
        matrix[n_sp,n_sp]Vphy;     // phylogeny
		
	}

transformed data{
	matrix[n_sp, n_sp] Ident; 
        Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
}

parameters {
  real<lower=0> vsigma; 
  real<lower=0,upper=10> lam;  
  real a; // intercept
  real b_force; // slope of forcing effect 
	}

model {
       vector[N] yhat;
       matrix[n_sp, n_sp] vlambda;
       matrix[n_sp, n_sp] sigma_y;   
       vlambda = (lam*Lmat + Ident) .* Vphy;
       sigma_y = vlambda*vsigma; 

       yhat = a + b_force * x;
        
	a ~ normal(0,10);
        b_force ~ normal(0,3);
	lam ~ normal(0, 5);
	vsigma ~ normal(0, 10);
	y ~ multi_normal(yhat, sigma_y);

}

