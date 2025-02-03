// Feb 3, 2025
// Started by V Van der Meersch
// Inspired by previous works on PPMs in the lab

// Beta regression with proportion and degenerate responses (i.e. 0 and 1)
// Original code from R. Kubinec, https://github.com/saudiwin/ordbetareg

functions {
  
  // prior from Michael Betancourt for ordered cutpoints
    real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data {
  
  int<lower=0> N_prop; // number of proportion observations (0,1)
  int<lower=0> N_degen; // number of 0/1 observations
  
  int<lower=1> Nsp; // number of species
  int<lower=1, upper=Nsp> sp_prop[N_prop]; // species ID for prop. observations
  int<lower=1, upper=Nsp> sp_degen[N_degen]; // species ID for degen. observations
  
  vector[N_prop] y_prop; // Y in (0,1)
  int<lower=0,upper = 1> y_degen[N_degen]; // Y in {0,1}
  
  vector[N_prop] x_prop; // covariate X for proportion outcome
  vector[N_degen] x_degen; // covariate X for degenerate (0,1) outcome
  
  corr_matrix[Nsp] Vphy; // phylogenetic relationship matrix (fixed)
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
  
  ordered[2] cutpoints; // cutpoints on ordered (latent) variable (also stand in as intercepts)
  real<lower=0> kappa; // scale parameter for beta regression
}

transformed parameters {
  
  corr_matrix[Nsp] C_a = lambda_a * Vphy;
  C_a = C_a - diag_matrix(diagonal(C_a)) + diag_matrix(diagonal(Vphy));
  
  corr_matrix[Nsp] C_b = lambda_b * Vphy;
  C_b = C_b - diag_matrix(diagonal(C_b)) + diag_matrix(diagonal(Vphy));
  
  // more numerically stable and more efficient to use pre-factored covariance matrices (i.e. multi_normal_cholesky in the following
  matrix[Nsp,Nsp] L_a = cholesky_decompose(sigma_a^2*C_a);
  matrix[Nsp,Nsp] L_b =  cholesky_decompose(sigma_b^2*C_b); 
  
  real calc_degen[N_degen];
  real calc_prop[N_prop];
  
  if(N_degen>0) {
    for(i in 1:N_degen){
      calc_degen[i] = a[sp_degen[i]] + b[sp_degen[i]] * x_degen[i];
    }
  }
  
  for(i in 1:N_prop){
    
    calc_prop[i] = a[sp_prop[i]] + b[sp_prop[i]] * x_prop[i];
    
  }

}
model {
  
  a ~ multi_normal_cholesky(rep_vector(a_z,Nsp), L_a); 
  b ~ multi_normal_cholesky(rep_vector(b_z,Nsp), L_b); 
  
  target += induced_dirichlet_lpdf(cutpoints | rep_vector(1, 3), 0);
  
  // need separate loops for logit (0/1) and beta regression
  if(N_degen>0) {
    for(n in 1:N_degen) {
      if(y_degen[n]==0) {
        // Pr(Y==0)
        target += log1m_inv_logit(calc_degen[n] - cutpoints[1]);
      } else {
        //Pr(Y==1)
        target += log_inv_logit(calc_degen[n] - cutpoints[2]);
      }
    }
  }
  
  for(n in 1:N_prop) {
    // Pr(Y in (0,1))
    target += log(inv_logit(calc_prop[n] - cutpoints[1]) - inv_logit(calc_prop[n] - cutpoints[2]));
    // Pr(Y==x where x in (0,1))
    y_prop[n] ~ beta_proportion(inv_logit(calc_prop[n]),kappa);
  }
  
  // priors
  a_z ~ normal(0, 1.5); 
  b_z ~ normal(0.5, 1); 
  
  lambda_a ~ beta(1.5, 1.5);
  sigma_a ~ normal(0, 1);
  
  lambda_b ~ beta(1.5, 1.5);
  sigma_b ~ normal(0, 1);
  
  kappa ~ exponential(.1); // vague?
  
}


