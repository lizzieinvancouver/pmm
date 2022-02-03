## Started 13 Jan 2021 ##
## Adapted from earlier code (ubermini_2.R)
## Generate test data for phylogeny Stan model ##

## Modified by Deirdre (Dec 2021) to modify the code to increase model efficiency for large matricies.

## Feb 2: I think we are close - want to double check that this model, with Geoff's edits and no phylogeny on the intercept
# is actually giving us good values

if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("Lizzie", getwd())>0)) {
    setwd("~/Documents/git/teaching/stan/pmm/analyses")
} else{
    setwd("/home/deirdre/pmm") # for midge
}

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)
require(tidybayes)

# Set number of cores available
## options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

# Set seed
set.seed(2021)

nspecies = 800
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
param <- list(
     mu_a = 4, # root value intercept
     sigma_a = 1,
    #           lam_interceptsa = 0.4, # lambda intercept
    #           sigma_interceptsa = 0.2, # rate of evolution intercept
              b_zf = 0.6, # root value trait1 slope
              lam_interceptsbf = 0.7, # lambda trait1
              sigma_interceptsbf = 0.1, # rate of evolution trait1
              sigma_y = 0.01 # overall sigma
              )
# Set priors
phypriors <- list(
    mu_a_prior_mu = 4, # true value
    mu_a_prior_sigma = 1,
    sigma_a_prior_mu = 0.5,
    sigma_a_prior_sigma = 1,
    # lam_interceptsa_prior_alpha = 4, # 
    # lam_interceptsa_prior_beta = 6, # 
    # sigma_interceptsa_prior_mu = 0.2, # true value
    # sigma_interceptsa_prior_sigma = 0.2,
    b_zf_prior_mu = 0.6, # true value
    b_zf_prior_sigma = 1,
    lam_interceptsbf_prior_alpha = 7, #
    lam_interceptsbf_prior_beta = 3, # 
    sigma_interceptsbf_prior_mu = 0.1, # true value
    sigma_interceptsbf_prior_sigma = 0.1,
    sigma_y_prior_mu = 0.01,
    sigma_y_prior_sigma = 1
    
   
)

# Generate intercept
intercepts <- rnorm(nspecies, 10, param[["sigma_a"]])

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
    # a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = 1)
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_zf"]], sd = 1)
    return(append(list(b_force = b_z.temp),
                  param))
}

# Fit model

testme <- stan("stan/uber_oneSlope_noIntercept_cholesky.stan",
               data = append(list(N=nrow(dfhere),
                                  n_sp=nspecies,
                                  sp=dfhere$sp,
                                  x1=dfhere$x1,
                                  y=dfhere$y,
                                  Vphy=vcv(spetree, corr = TRUE)),
                             phypriors),
               init = simu_inits,
               iter = 2000,
               warmup = 1000,
               chains = 4,
               seed = 62921
               )
save(testme, file = "output_phylo_cholesky_slopenoint.Rda")

summary(testme)$summary[c("b_zf","lam_interceptsbf","sigma_interceptsbf","sigma_y"),"mean"]; t(param)
# 


# 
# load("output/phylo_cholesky_benmod_Jan25.Rda")
# # 
#  ssm<- as.shinystan(test_new)
#  launch_shinystan(ssm)

# get_variables(mdl.b)
# summary(mdl.b)$summary[c( "b_zf","lam_interceptsbf","sigma_interceptsbf","sigma_y"),"mean"]; t(param)

# pairs(test_new, pars = c("a_z","lam_interceptsa","sigma_interceptsa", "b_zf","lam_interceptsbf","sigma_interceptsbf","sigma_y", "lp__")) 
# Summarize fit
# trying to get the ESS based on code from: 
# https://betanalpha.github.io/assets/case_studies/identifiability.html
# step times plot:
# stepsizes1 <- sapply(1:4, function(c) get_sampler_params(test_new, inc_warmup=FALSE)[[c]][,'stepsize__'][1])
# steps1 <- do.call(rbind, get_sampler_params(test_new, inc_warmup=FALSE))[,'n_leapfrog__']
# table(steps1)
# int_times1 <- unlist(lapply(1:4, function(c) stepsizes1[c] * steps1[(1000 * (c - 1) + 1): (1000 * c)]))
# 
# int_time_breaks <- seq(-0.1, 8.1, 0.1)
# 
# hist(int_times1, breaks=int_time_breaks, main="N = 1",
#      col="darkgreen", border="blue",
#      xlab="Integration Time", yaxt='n', ylab="")
# 
# summary1 <- summary(test_new, probs = c(0.5))$summary
# 
# sigma_y_ess <- summary1[,'n_eff'][1]
# lam_inta_ess <- summary1[,'n_eff'][2]
# sig_inta_ess <- summary1[,'n_eff'][3]
# lam_intb_ess <- summary1[,'n_eff'][4]
# sig_intb_ess <- summary1[,'n_eff'][5]
# b_z_ess <- summary1[,'n_eff'][6]
# a_z_ess <- summary1[,'n_eff'][7]
# 
# ess <- summary1[1:6,'n_eff']
# 
# 
# # ESS per leapfrog step - measure of computational efficiency, decreases w/ no. obs, as cost of leapfrogs incre ESS per computational cost decreases
# ess_per <- c(summary1[1:6,'n_eff'] / sum(steps1))
# #b_ess_per <- c(summary1[,'n_eff'][2] / sum(steps1))
# 
# par(mfrow=c(1, 1))
# 
# plot(1, a_ess_per / a_ess_per[1], col="darkgreen", type="l", log="y",
#      xlab="N", xaxt='n',
#      ylab="Relative Effective Sample Size / Grad Eval (a)", ylim=c(0.01, 1))
# axis(1, at=1:3, labels=c("1", "100", "10000"))
# 
# plot(1, b_ess_per / b_ess_per[1], col="darkgreen", type="l", log="y",
#      xlab="N", xaxt='n',
#      ylab="Relative Effective Sample Size / Grad Eval (b)", ylim=c(0.01, 1))
# axis(1, at=1:3, labels=c("1", "100", "10000"))

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
