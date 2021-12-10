## Started 13 Jan 2021 ##
## Adapted from earlier code (ubermini_2.R)
## Generate test data for phylogeny Stan model ##

## Modified by Deirdre (Dec 2021) to modify the code to increase model efficiency for large matricies.

if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("Lizzie", getwd())>0)) {
    setwd("~/Documents/git/teaching/stan/pmm/analyses")
} else{
    setwd("/home/deirdre/phylogeny") # for midge
}

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

nspecies = 200
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
param <- list(a_z = 4, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              b_zf = 0.6, # root value trait1 slope
              lam_interceptsbf = 0.7, # lambda trait1
              sigma_interceptsbf = 0.1, # rate of evolution trait1
              sigma_y = 0.01 # overall sigma
              )
# Set priors
phypriors <- list(
    a_z_prior_mu = 4, # true value
    a_z_prior_sigma = 1,
    lam_interceptsa_prior_alpha = 4, # 
    lam_interceptsa_prior_beta = 6, # 
    sigma_interceptsa_prior_mu = 0.2, # true value
    sigma_interceptsa_prior_sigma = 0.2,
    b_zf_prior_mu = 0.6, # true value
    b_zf_prior_sigma = 1,
    lam_interceptsbf_prior_alpha = 7, #
    lam_interceptsbf_prior_beta = 3, # 
    sigma_interceptsbf_prior_mu = 0.1, # true value
    sigma_interceptsbf_prior_sigma = 0.1,
    sigma_y_mu_prior = 0.01, # true value
    sigma_y_mu_sigma = 0.01,
    mu_prior_a_ = 0, # adding priors for a_ and b_
    sigma_prior_a_ = 1,
    mu_prior_b_ = 0, 
    sigma_prior_b_ = 1
    
)

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsbf"]])         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)

dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       intercept = rep(intercepts[i], nind),
                       x1 = rnorm(n = nind, mean = 10, sd = 3),
                       trait1 = rep(slopes_bf[i], nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$intercept + dfhere$x1 * dfhere$trait1
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = 1)
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_zf"]], sd = 1)
    return(append(list(a = a_z.temp, b_force = b_z.temp),
                  param))
}

# Fit model

# test_old <- stan("Stan/uber_oneslopeintercept_cholesky.stan",
#                data = append(list(N=nrow(dfhere),
#                                   n_sp=nspecies,
#                                   sp=dfhere$sp,
#                                   x1=dfhere$x1,
#                                   y=dfhere$y,
#                                   Vphy=vcv(spetree, corr = TRUE)),
#                              phypriors),
#                init = simu_inits,
#                iter = 4000,
#                warmup = 2000,
#                chains = 4,
#                seed = 62921
#                )
# 
# save(test_old, file = "output_phylo_cholesky_oldchol.Rda")

test_new <- stan("analyses/stan/uber_oneslopeintercept_cholesky_modified.stan",
                 data = append(list(N=nrow(dfhere),
                                    n_sp=nspecies,
                                    sp=dfhere$sp,
                                    x1=dfhere$x1,
                                    y=dfhere$y,
                                    Vphy=vcv(spetree, corr = TRUE)),
                               phypriors),
                 #init = simu_inits,
                 iter = 4000,
                 warmup = 3000,
                 chains = 4)

save(test_new, file = "cholesky_fixed_noint.Rda")

#load("analyses/output/cholesky_fixed_noint.Rda")
load("analyses/output/cholesky_fixed.Rda")

ssm<- as.shinystan(test_new)
launch_shinystan(ssm)

# Summarize fit
summary(test_new)$summary[c("a_z","lam_interceptsa","sigma_interceptsa", "b_zf","lam_interceptsbf","sigma_interceptsbf","sigma_y"),"mean"]; t(param)
# 
# # Compare to true values
# summary(testme)$summary[names(param), "mean"]
# t(param)
# ## summary(testme)$summary[names(param), "2.5%"]
# ## summary(testme)$summary[names(param), "97.5%"]
# ## summary(testme)$summary[names(param), "50%"]
# 
# # Compare to GEIGER-fitted model
# fitContinuous(spetree, intercepts, model="lambda")
# fitContinuous(spetree, slopes_bf, model="lambda")
