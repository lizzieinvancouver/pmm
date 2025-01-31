// Jan 24, 2025
// Started by V Van der Meersch
// Inspired by previous works on PPMs in the lab

data {
  
  int<lower=1> N;
  int<lower=1> Nsp;
  int<lower=1, upper=Nsp> sp[N];
  vector[N] y; // response
  vector[N] x; // predictor
  corr_matrix[Nsp] Vphy; // phylogenetic relationship matrix
  
}

parameters {
  
  // intercept (+ the portion of phenotypes, not predicted by x, which still covary between related species, right?)
  vector[Nsp] a; 
  real a_z; // root value
  real<lower=0, upper=1> lambda_a; // phylogenetic structure      
  real<lower=0> sigma_a; // overall rate of change (brownian motion?)
  
  // slope of x effect
  vector[Nsp] b; 
  real b_z; // root value
  real<lower=0, upper=1> lambda_b;  // phylogenetic structure        
  real<lower=0> sigma_b; // overall rate of change (brownian motion?)
  
  real<lower=0> sigma_y; // others sources of error
  
}

transformed parameters{
  
  corr_matrix[Nsp] C_a = lambda_a * Vphy;
  C_a = C_a - diag_matrix(diagonal(C_a)) + diag_matrix(diagonal(Vphy));
  
  corr_matrix[Nsp] C_b = lambda_b * Vphy;
  C_b = C_b - diag_matrix(diagonal(C_b)) + diag_matrix(diagonal(Vphy));
  
  // more numerically stable and more efficient to use pre-factored covariance matrices (i.e. multi_normal_cholesky in the following
  matrix[Nsp,Nsp] L_a = cholesky_decompose(sigma_a^2*C_a); 
  matrix[Nsp,Nsp] L_b =  cholesky_decompose(sigma_b^2*C_b); 
  
}

model {
  
  real yhat[N];
  
  for(i in 1:N){
    yhat[i] = a[sp[i]] + b[sp[i]] * x[i];
  }
  
  a ~ multi_normal_cholesky(rep_vector(a_z,Nsp), L_a); 
  b ~ multi_normal_cholesky(rep_vector(b_z,Nsp), L_b); 
  
  y ~ normal(yhat, sigma_y);
  
  
  
  // priors (from phyloMdlLambdaIntSlope.stan)
  a_z ~ normal(30, 10); 
  b_z ~ normal(-2, 10); 

  lambda_a ~ beta(1, 1);
  sigma_a ~ normal(30, 20);
  
  lambda_b ~ beta(1, 1);
  sigma_b ~ normal(1, 5);
    
  sigma_y ~ normal(10, 10); 
  
}

