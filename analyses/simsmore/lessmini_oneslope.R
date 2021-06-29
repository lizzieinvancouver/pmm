## Started 13 Jan 2021 ##
## Adapted from earlier code (ubermini_2.R)
## Generate test data for phylogeny Stan model ##

## load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)
library(truncnorm)

## options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

nspecies = 90
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait parameters
param <- list(b_z = 0.6, # root value trait1
              lam_interceptsb = 0.7, # lambda trait1
              sigma_interceptsb = 0.1, # rate of evolution trait1
              sigma_y = 0.01 # overall sigma
              )

# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsb"]])         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_z"]], mu = 0, sig2 = param[["sigma_interceptsb"]] ^ 2)

dfhere <- data.frame(sp = c(), x1=numeric(), trait1=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       x1 = rnorm(n = nind, mean = 10, sd = 3),
                       trait1 = rep(slopes_bf[i], nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$x1 * dfhere$trait1
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_z"]], sd = 1)
    return(append(list(b_force = b_z.temp),
                  param))
}

# Fit model
testme <- stan("lessmini_oneslope.stan",
               data=list(N=nrow(dfhere),
                         n_sp=nspecies,
                         sp=dfhere$sp,
                         x1=dfhere$x1,
                         y=dfhere$y,
                         Vphy=vcv(spetree, corr = TRUE)),
               iter=4000, chains=4,
               init = simu_inits
               ) # seed=123456

# Summarize fit
summary(testme)$summary

# Compare to true values
summary(testme)$summary[names(param), "mean"]
t(param)
## summary(testme)$summary[names(param), "2.5%"]
## summary(testme)$summary[names(param), "97.5%"]
## summary(testme)$summary[names(param), "50%"]

fitContinuous(spetree, slopes_bf, model="lambda")
