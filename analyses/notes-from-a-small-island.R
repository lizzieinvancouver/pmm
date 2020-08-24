# Headers
library(rstan)
library(geiger)
library(phytools)

# Classic PGLS-type estimation in stan (?)
signal.stan <- "data
{
   int n_spp;
   matrix[n_spp,n_spp] Vphy_lam;
   vector[n_spp] trait;
   vector[n_spp] diag;
}
parameters
{
  real<lower=0> sigma;
  real<lower=0> lam;
  real z;
}
model
{
  sigma ~ normal(0,1);
  lam ~ normal(0,1);
  z ~ normal(0, .5);
  trait ~ multi_normal(
    rep_vector(z,n_spp), (diag_matrix(diag) + lam*Vphy_lam) * sigma
  );
}
"
signal.model <- stan_model(model_code=signal.stan)

# Simulate some data and then check
signal.sim <- function(lam, n_spp=100){
    data <- list(n_spp=n_spp)
    data$tree <- sim.bdtree(n=data$n_spp)
    data$Vphy_lam <- data$Vphy <- vcv(data$tree)
    data$diag <- diag(data$Vphy)
    diag(data$Vphy_lam) <- 0
    data$trait <- fastBM(rescale(data$tree, "lambda", lam))
    return(data)
}
data.zero <- signal.sim(0)
data.half <- signal.sim(.5)
data.one <- signal.sim(1)

.wrap.sig <- function(x, ...){
    model <- sampling(signal.model, x, ...)
    print(setNames(
        c(summary(model.zero)$summary["lam","mean"], fitContinuous(data.zero$tree, data.zero$trait, model="lambda")$opt$lambda),
        c("stan","geiger")
    ))
    hist(unlist(extract(model.zero, "lam")))
    return(model)
}
model.zero <- .wrap.sig(data.zero, iter=3000, warmup=2500, chains=2)
model.half <- .wrap.sig(data.half, iter=3000, warmup=2500, chains=2)
model.one <- .wrap.sig(data.one, iter=3000, warmup=2500, chains=2)
# ... hope you like divergent transitions! ...
