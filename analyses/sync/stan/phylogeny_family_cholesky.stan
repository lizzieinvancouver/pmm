// Aloha!
// By Megan & Lizzie started on 6 January 2017
// Modified by Deirdre for her synchrony analysis - adding family level phylogenetic effects
// Trying to do a simple 3-level model
// Imagine some data: 
    // You have observed densities of something across plots, and plots occur within sites
    // DL: or in my case you have the nested effect of species within a family
    // 3 level: 
        // (1) Observation level (density) - or indiviudals day of a phenological event
        // (2) plot level (multiple observations taken within each plot) - or species
        // (3) Multiple plots for each site - or taxonomic family
    // Predictor data (aka x) observed only at observation level (imagine you measured moisture around each plant or such)
    // ... (that is, no additional predictors added at plot or site level)
// Model structure does not assume same number of plots across sites
// Code below includes pooled intercepts and slopes at all levels
// A link (we did not use it, but I am keeping it here for now): http://stackoverflow.com/questions/29379001/nested-model-in-stan
functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
    return(quad_form_diag(local_vcv, rep_vector(sigma, rows(vcv))));
  }
}

data {
  int N;    //total number of observations (e.g., N=120=3 observations in each plot x 4 plots per site x 10 sites )
  int Nspp;    //total number of plots/species 
  int Nfam;    //total number of sites/families
  vector[N] y;    // the data, y, is a vector of length N
  int sppnum[N];  // column of plot number identifiers for each obs
  int famnum[Nspp];  // column of site number identifiers for each plot
  vector[N] x;        // vector of predictor X for each obs
  matrix[Nfam,Nfam] Vphy;     // phylogeny
  
  //priors
  real b_z_prior_mu;
  real b_z_prior_sigma;
  real lam_interceptsb_prior_alpha;
  real lam_interceptsb_prior_beta;
  real sigma_interceptsb_prior_mu;
  real sigma_interceptsb_prior_sigma;
  real sigma_y_mu_prior;
  real sigma_y_sigma_prior;  
}

parameters {
  vector[Nspp] a_spp;    
  vector[Nspp] b_spp;    
  vector[Nfam] a_fam;    
  vector[Nfam] b_fam;    
  real b_z;
  
  real<lower=0> sig_a_fam;  //variance in intercept across families; 
  real<lower=0> sig_b_fam;  //variance in slopes across families; 
  real mu_a;                    //mean intercept across species; 
      // the species intercept for species s is drawn from distribution with mean mu_a...
  real<lower=0> sig_a;          //...and standard deviation sig_a
  real mu_b;                    //mean slope across species; 
      // the species slope are drawn from distribution with mean mu_b...
  //real<lower=0> sig_b;          //...and standard deviation sig_b
  real<lower=0> sig_y;         // observation error
  
  real<lower=0, upper=1> lam_interceptsb;       
  real<lower=0> sigma_interceptsb;    
}

model {
  real yhat[N];                 //vector of predicted y's (`ypred' for Lizzie)
  
  matrix[Nfam, Nfam] vcv_b; 
  
  for (i in 1:N) {
    yhat[i] = a_spp[sppnum[i]] + b_spp[sppnum[i]]*x[i];
  }

  //For estimating a single value for all within-site variacnes
  for (j in 1:Nspp){
    a_spp[j] ~ normal(a_fam[famnum[j]], sig_a_fam);
    b_spp[j] ~ normal(b_fam[famnum[j]], sig_b_fam);
  }
  
  y ~ normal(yhat,sig_y);        //data is distributed normally around predicted (yhat) with s.d. sig_y (this is error of data around predicted values)
    
  vcv_b = cholesky_decompose(lambda_vcv(Vphy, lam_interceptsb, sigma_interceptsb));
  
  a_fam ~ normal(mu_a,sig_a);
  b_fam ~ multi_normal_cholesky(rep_vector(b_z, Nfam), vcv_b);
  
  sig_a ~ normal(0.5,1);
  mu_a ~ normal(0, sig_a);
 
  b_z ~ normal(b_z_prior_mu, b_z_prior_sigma);
  
  lam_interceptsb ~ beta(lam_interceptsb_prior_alpha, lam_interceptsb_prior_beta);
  sigma_interceptsb ~ normal(sigma_interceptsb_prior_mu, sigma_interceptsb_prior_sigma);
  sig_y ~ normal(sigma_y_mu_prior, sigma_y_sigma_prior);
  

  

}

