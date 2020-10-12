functions {
  matrix vcv_mat(matrix corr_mat, real sigma){
    matrix[rows(corr_mat),cols(corr_mat)] phylo_cov;
    matrix[rows(corr_mat),cols(corr_mat)] sigma_mat;
    sigma_mat = diag_matrix(rep_vector(sigma, rows(corr_mat)));
    phylo_cov = sigma_mat * corr_mat * sigma_mat;
    return(phylo_cov);
  }
}

data {
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] x; 	// predictor
        matrix[n_sp,n_sp]Vphy;     // phylogeny
}

parameters {
  real b_force_latent_mu; // mean forcing slope before phylogenetic effect
  vector[n_sp] b_force_latent;
  vector[n_sp] b_force_phylogeny;
  real<lower=0> sigma_y;
  real<lower=0> sigma_phylogeny;
  real<lower=0> b_force_latent_sigma; // variation in slope before phylogenetic effect
}

transformed parameters {
  real yhat[N];
 
  for(i in 1:N){
    yhat[i] =  x[i] * (b_force_latent[sp[i]] + b_force_phylogeny[sp[i]]);
  }
}

model {
  y ~ normal(yhat, sigma_y);
  
  b_force_latent ~ normal(b_force_latent_mu, b_force_latent_sigma);
  b_force_phylogeny ~ multi_normal(rep_vector(0, n_sp), vcv_mat(Vphy, sigma_phylogeny));
  
  b_force_latent_mu ~ normal(2, 2);
  b_force_latent_sigma ~ normal(0, 2);
  sigma_phylogeny ~ normal(0, 1);
  sigma_y ~ normal(0, 2);
  
}

