## Adapted from earlier code (ubermini_2.R)
## Generate test data for phylogeny Stan model ##

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

# Set number of cores available
## options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

# Set seed
set.seed(2021)

nspecies = 35
nind = 5

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
param <- list(a_z = 4, # root value intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              b_zf = 0.6, # root value trait1
              sigma_interceptsbf = 0.1, # rate of evolution trait1
              b_zc = 0.88, # root value trait2
              sigma_interceptsbc = 0.07, # rate of evolution trait2
              b_zp = 1.1, # root value trait3
              sigma_interceptsbp = 0.3, # rate of evolution trait3
              sigma_y = 0.01 # overall sigma
              )
# Set priors
phypriors <- list(
    a_z_prior_mu = 4, # true value
    a_z_prior_sigma = 1,
    sigma_interceptsa_prior_mu = 0.2, # true value
    sigma_interceptsa_prior_sigma = 0.2,
    b_zf_prior_mu = 0.6, # true value
    b_zf_prior_sigma = 1,
    sigma_interceptsbf_prior_mu = 0.1, # true value
    sigma_interceptsbf_prior_sigma = 0.1,
    b_zc_prior_mu = 0.6, # true value
    b_zc_prior_sigma = 1,
    sigma_interceptsbc_prior_mu = 0.07, # true value
    sigma_interceptsbc_prior_sigma = 0.07,
    b_zp_prior_mu = 1.1, # true value
    b_zp_prior_sigma = 1,
    sigma_interceptsbp_prior_mu = 0.3, # true value
    sigma_interceptsbp_prior_sigma = 0.1,
    sigma_y_mu_prior = 0.01, # true value
    sigma_y_mu_sigma = 0.01
)

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", 0)         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", 0)         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)
# Generate bc slope
scaledtree_bc <- rescale(spetree, model = "lambda", 0)         
slopes_bc <- fastBM(scaledtree_bc, a = param[["b_zc"]], mu = 0, sig2 = param[["sigma_interceptsbc"]] ^ 2)
# Generate bp slope
scaledtree_bp <- rescale(spetree, model = "lambda", 0)         
slopes_bp <- fastBM(scaledtree_bp, a = param[["b_zp"]], mu = 0, sig2 = param[["sigma_interceptsbp"]] ^ 2)

dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric(), x2 = numeric(), trait2 = numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       intercept = rep(intercepts[i], nind),
                       x1 = rnorm(n = nind, mean = 0, sd = 5),
                       trait1 = rep(slopes_bf[i], nind),
                       x2 = rnorm(n = nind, mean = 0, sd = 5),
                       trait2 = rep(slopes_bc[i], nind),
                       x3 = rnorm(n = nind, mean = 0, sd = 5),
                       trait3 = rep(slopes_bp[i], nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$intercept + dfhere$x1 * dfhere$trait1 + dfhere$x2 * dfhere$trait2 + dfhere$x3 * dfhere$trait3
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = .5)
    b_zf.temp <- rnorm(n = nspecies, mean = param[["b_zf"]], sd = .5)
    b_zc.temp <- rnorm(n = nspecies, mean = param[["b_zc"]], sd = .5)
    b_zp.temp <- rnorm(n = nspecies, mean = param[["b_zp"]], sd = .5)
    return(append(list(a = a_z.temp,
                       b_force = b_zf.temp,
                       b_chill = b_zc.temp,
                       b_photo = b_zp.temp),
                  param))
}

# Fit model
testme <- stan("stan/uber_threeslopeintercept_lambda0.stan",
               data = append(list(N=nrow(dfhere),
                                  n_sp=nspecies,
                                  sp=dfhere$sp,
                                  x1=dfhere$x1,
                                  x2=dfhere$x2,
                                  x3=dfhere$x3,
                                  y=dfhere$y,
                                  Vphy=vcv(spetree, corr = TRUE)),
                             phypriors),
               init = simu_inits,
               iter = 2000,
               warmup = 1000,
               chains = 4,
               seed = 62921,
               refresh = 20
               )

# Summarize fit
summary(testme)$summary

# Compare to true values
summary(testme)$summary[names(param), "mean"]
t(param)

# Compare to GEIGER-fitted model
fitContinuous(spetree, intercepts, model="lambda")
fitContinuous(spetree, slopes_bf, model="lambda")
fitContinuous(spetree, slopes_bc, model="lambda")
fitContinuous(spetree, slopes_bp, model="lambda")
