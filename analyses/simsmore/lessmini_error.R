## Started 9 July 2021 ##
## Adapted from lessmini_oneslope.r ##
## Take the phylo structure off slope and put on extra 'error' term ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory.
# Add in your own path in an if statement for your file structure (please don't delete this)
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/teaching/stan/pmm/analyses/simsmore") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("heregoboomboom")

## load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

## options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

nspecies = 90
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait parameters
param <- list(b_z = 0.6, # root value trait1 (leave as is, but think of it as this extra 'error')
              lam_interceptsb = 0.7, # lambda trait1
              sigma_interceptsb = 0.1, # rate of evolution trait1
              sigma_y = 0.01, # overall sigma
              slope_mu = 3, # new slope mean
              slope_sigma = 1 # new slope variance
              )

# Set priors
phypriors <- list(
    b_z_prior_mu = 1,
    b_z_prior_sigma = 3,
    lam_interceptsb_prior_alpha = 7,
    lam_interceptsb_prior_beta = 3,
    sigma_interceptsb_prior_mu = 0.1,
    sigma_interceptsb_prior_sigma = 1,
    sigma_y_mu_prior = 1, 
    sigma_y_mu_sigma = 5
)

# Generate slopes -- no phylo structure
slopez <- rnorm(nspecies, param[["slope_mu"]], param[["slope_sigma"]])

# Generate bf 
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsb"]])         
error_bf <- fastBM(scaledtree_bf, a = param[["b_z"]], mu = 0, sig2 = param[["sigma_interceptsb"]] ^ 2)

dfhere <- data.frame(sp = c(), x1=numeric(), trait1=numeric(), slope=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       x1 = rnorm(n = nind, mean = 10, sd = 3),
                       trait1 = rep(error_bf[i], nind),
                       slope = rep(slopez[i], each=nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$x1 * dfhere$slope + dfhere$trait1 # adjusted here
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_z"]], sd = 1)
    return(append(list(b_force = b_z.temp),
                  param))
}

# Fit model
testme <- stan("..//stan/uber_oneslope.stan",
               data=append(list(N=nrow(dfhere),
                         n_sp=nspecies,
                         sp=dfhere$sp,
                         x1=dfhere$x1,
                         y=dfhere$y,
                         Vphy=vcv(spetree, corr = TRUE)),
               phypriors),
               init = simu_inits,
               iter = 3000,
               warmup = 2000,
               chains = 4)

# Summarize fit
summary(testme)$summary

# Compare to true values
summary(testme)$summary[names(param)[1:4], "mean"]
t(param)
## summary(testme)$summary[names(param), "2.5%"]
## summary(testme)$summary[names(param), "97.5%"]
## summary(testme)$summary[names(param), "50%"]

fitContinuous(spetree, error_bf, model="lambda")
