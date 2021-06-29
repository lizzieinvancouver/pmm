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
nind = 20

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait parameters
param <- list(a_z = 4, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              b_zf = 0.6, # root value trait1
              lam_interceptsbf = 0.7, # lambda trait1
              sigma_interceptsbf = 0.1, # rate of evolution trait1
              b_zc = 0.88, # root value trait2
              lam_interceptsbc = 0.9, # lambda trait2
              sigma_interceptsbc = 0.07, # rate of evolution trait2
              b_zp = 1.1, # root value trait3
              lam_interceptsbp = 0.45, # lambda trait3
              sigma_interceptsbp = 0.3, # rate of evolution trait3
              sigma_y = 0.02 # overall sigma
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
testme <- stan("lessmini_threeslopeintercept.stan",
               data=list(N=nrow(dfhere),
                         n_sp=nspecies,
                         sp=dfhere$sp,
                         x1=dfhere$x1,
                         x2=dfhere$x2,
                         x3=dfhere$x3,
                         y=dfhere$y,
                         Vphy=vcv(spetree, corr = TRUE)),
               iter=2000, chains=4,
               init = simu_inits
               ) # seed=123456

# Summarize fit
head(summary(testme)$summary)

# Compare to true values
summary(testme)$summary[names(param), "mean"]
t(param)
## summary(testme)$summary[names(param), "2.5%"]
## summary(testme)$summary[names(param), "97.5%"]
## summary(testme)$summary[names(param), "50%"]

fitContinuous(spetree, slopes_bf, model="lambda")
fitContinuous(spetree, intercepts, model="lambda")

## plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])
## abline(a = 0, b = 1, lty = "dotted")
