## Started 13 Jan 2021 ##
## Adapted from earlier code (ubermini_2.R)
## Generate test data for phylogeny Stan model ##

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

setwd("~/Documents/git/teaching/stan/pmm/analyses/")

# Set number of cores available
## options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

# Set seed
set.seed(2021)

nspecies = 40
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
param <- list(a_z = 30, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 16, # rate of evolution intercept
              b_zf = -6, # root value trait1
              lam_interceptsbf = 0.7, # lambda trait1
              sigma_interceptsbf = 6, # rate of evolution trait1
              b_zc = -7, # root value trait2
              lam_interceptsbc = 0.5, # lambda trait2
              sigma_interceptsbc = 7, # rate of evolution trait2
              b_zp = -1.4, # root value trait3
              lam_interceptsbp = 0.3, # lambda trait3
              sigma_interceptsbp = 2, # rate of evolution trait3
              sigma_y = 13 # overall sigma
              )
# Set priors
phypriors <- list(
    a_z_prior_mu = 30, 
    a_z_prior_sigma = 10,
    lam_interceptsa_prior_alpha = 6, 
    lam_interceptsa_prior_beta = 6,  
    sigma_interceptsa_prior_mu = 10, 
    sigma_interceptsa_prior_sigma = 10,
    b_zf_prior_mu = 0, 
    b_zf_prior_sigma = 20, 
    lam_interceptsbf_prior_alpha = 6, 
    lam_interceptsbf_prior_beta = 6,  
    sigma_interceptsbf_prior_mu = 10, 
    sigma_interceptsbf_prior_sigma = 10,
    b_zc_prior_mu = 0, 
    b_zc_prior_sigma = 20,
    lam_interceptsbc_prior_alpha = 6, 
    lam_interceptsbc_prior_beta = 6,  
    sigma_interceptsbc_prior_mu = 10,
    sigma_interceptsbc_prior_sigma = 10,
    b_zp_prior_mu = 0,
    b_zp_prior_sigma = 20,
    lam_interceptsbp_prior_alpha = 6, 
    lam_interceptsbp_prior_beta = 6,  
    sigma_interceptsbp_prior_mu = 10, 
    sigma_interceptsbp_prior_sigma = 10,
    sigma_y_mu_prior = 10, 
    sigma_y_mu_sigma = 10
)

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsbf"]])         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)
# Generate bc slope
scaledtree_bc <- rescale(spetree, model = "lambda", param[["lam_interceptsbc"]])         
slopes_bc <- fastBM(scaledtree_bc, a = param[["b_zc"]], mu = 0, sig2 = param[["sigma_interceptsbc"]] ^ 2)
# Generate bp slope
scaledtree_bp <- rescale(spetree, model = "lambda", param[["lam_interceptsbp"]])         
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
testme <- stan("stan/uber_threeslopeintercept.stan",
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
# summary(testme)$summary

# Compare to true values
summary(testme)$summary[names(param), c("2.5%", "25%", "mean", "75%", "97.5%")]
t(param)
## summary(testme)$summary[names(param), "2.5%"]
## summary(testme)$summary[names(param), "97.5%"]
## summary(testme)$summary[names(param), "50%"]

# Compare to GEIGER-fitted model
fitContinuous(spetree, intercepts, model="lambda")
fitContinuous(spetree, slopes_bf, model="lambda")
fitContinuous(spetree, slopes_bc, model="lambda")
fitContinuous(spetree, slopes_bp, model="lambda")
