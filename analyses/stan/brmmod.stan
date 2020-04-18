// generated with brms 2.6.0'
// make_stancode of brmmod in run_MCMCbrms_pglmm.R 

functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  vector[N] Y;  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1
  int<lower=1> J_1[N];
  int<lower=1> N_1;
  int<lower=1> M_1;
  // cholesky factor of known covariance matrix 
  matrix[N_1, N_1] Lcov_1; 
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc = K - 1; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
} 
parameters { 
  vector[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
  real<lower=0> sigma;  // residual SD 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects
} 
transformed parameters { 
  // group-level effects 
  vector[N_1] r_1_1 = sd_1[1] * (Lcov_1 * z_1[1]);
} 
model { 
  vector[N] mu = temp_Intercept + Xc * b;
  for (n in 1:N) { 
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
  } 
  // priors including all constants 
  target += normal_lpdf(b | 0, 10); 
  target += normal_lpdf(temp_Intercept | 0, 50); 
  target += student_t_lpdf(sigma | 3, 0, 20)
    - 1 * student_t_lccdf(0 | 3, 0, 20); 
  target += student_t_lpdf(sd_1 | 3, 0, 20)
    - 1 * student_t_lccdf(0 | 3, 0, 20); 
  target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants 
  if (!prior_only) { 
    target += normal_lpdf(Y | mu, sigma);
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept - dot_product(means_X, b); 
  // additionally draw samples from priors
  real prior_b = normal_rng(0,10);
  real prior_sigma = student_t_rng(3,0,20);
  real prior_sd_1 = student_t_rng(3,0,20);
  // use rejection sampling for truncated priors
  while (prior_sigma < 0) {
    prior_sigma = student_t_rng(3,0,20);
  }
  while (prior_sd_1 < 0) {
    prior_sd_1 = student_t_rng(3,0,20);
  }
} 
