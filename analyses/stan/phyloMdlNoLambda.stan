// Oct 4, 2023
// Started by D Loughnan
//Aim of this model --- to check whether model with lambda out performs model without
// This is the NO lambda model--- just one cue (x)
//


data {
  int<lower=1> N;
  vector[N] yobs; 		// response
  vector[N] x1; 	// predictor (year)
  
  int<lower=1> Nspp;
  int<lower=1, upper= Nspp > species[N];
}


parameters {
  //real mu_grand;
  real<lower=0> sigma_y;    
  
  real mu_a; // grand mean
  vector[Nspp] a_sp; // intercept
  real<lower=0> sigma_a;
  
  real mu_b;
  real<lower=0> sigma_b;   
  
  
  vector[Nspp] b; 
  
  	}

	
model {
  vector[N] yhat = a_sp[species] + (( b[species])).* x1;

  a_sp ~ normal(mu_a, sigma_a);
  mu_a ~ normal(0,50);
  sigma_a ~ normal(0, 50);
  
  b ~ normal(mu_b, sigma_b);
  mu_b ~ normal(0,50);
  sigma_b ~ normal(0,10);
  
  sigma_y ~ normal(10,10);

  yobs ~ normal(yhat, sigma_y);
  
 
}









