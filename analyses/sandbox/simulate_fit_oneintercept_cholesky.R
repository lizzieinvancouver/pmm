## Started 18 May 2022 ##
## Copying from simulate_fit_oneslopeintercept_cholesky_dl.R to test grand mean 
## Generate test data for phylogeny Stan model ##

if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pmm/analyses/sync")
} else if(length(grep("lizzie", getwd())>0)) {
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

nspecies = 50
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
param <- list(a_z = 4, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              b_z = 0.6, # root value trait1 slope
              lam_interceptsb = 0.7, # lambda trait1
              sigma_interceptsb = 0.1, # rate of evolution trait1
              sigma_y = 0.01 # overall sigma
              )
# Set priors ... some of these are alarmingly narrow ...
# they need to be bigger before you can say your test data worked IMHO, but be careful to set properly for cholesky
phypriors <- list(
    a_z_prior_mu = 4, # true value
    a_z_prior_sigma = 3,
    lam_interceptsa_prior_alpha = 4, # 
    lam_interceptsa_prior_beta = 6, # 
    sigma_interceptsa_prior_mu = 0.2, # true value
    sigma_interceptsa_prior_sigma = 0.2,
    b_z_prior_mu = 0.6, # true value
    b_z_prior_sigma = 1,
    lam_interceptsb_prior_alpha = 7, #
    lam_interceptsb_prior_beta = 3, # 
    sigma_interceptsb_prior_mu = 0.1, # true value
    sigma_interceptsb_prior_sigma = 0.1,
    sigma_y_prior_mu = 0.01,
    sigma_y_prior_sigma = 1
    
   
)

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate b slope
scaledtree_b <- rescale(spetree, model = "lambda", param[["lam_interceptsb"]])         
slopes_b <- fastBM(scaledtree_b, a = param[["b_z"]], mu = 0, sig2 = param[["sigma_interceptsb"]] ^ 2)

dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       intercept = rep(intercepts[i], nind),
                       x1 = rnorm(n = nind, mean = 10, sd = 3),
                       trait1 = rep(slopes_b[i], nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$intercept + dfhere$x1 * dfhere$trait1
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = 1)
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_z"]], sd = 1)
    return(append(list(a = a_z.temp, b_force = b_z.temp),
                  param))
}

# Fit model

test <- stan("sandbox/stan/oneslopeinterceptcholforsync_grandmean.stan", # oneslopeinterceptcholforsync_grandmeans.stan
               data = append(list(N = nrow(dfhere),
                                  Nspp = nspecies,
                                  species = dfhere$sp,
                                  year = dfhere$x1,
                                  ypred = dfhere$y,
                                  Vphy = vcv(spetree, corr = TRUE)),
                             phypriors),
               # init = simu_inits,
               iter = 3000,
               warmup = 2000,
               chains = 4
               # seed = 62921
             )


# 19 May 2022: Notes on grand mean by Lizzie
# I made two different stan files (noted above) that pull out the grand means. In running them you can get 10-2000ish divergent transitions (depending on the seed I think) but I got the feeling the one pulling out just the intercept generally ran better than pulling out both. They always got close to the given values though. 
# 15 May 2022: General notes by Lizzie
# When I started running this originally, I was getting 70.68 seconds for the above with or without init turned off. Which is not a speed-up beyond way back when with ubermini_2.stan (I checked) ... And we end up in transposed units from the choleskey (which is annoying)
# Might be good to some day double-check this is true speed-up (same seeds etc.).
# We should also be careful to get to into initializing near the data always.

summary(test)$summary[c("a_z","lam_interceptsa","sigma_interceptsa",
                        "b_z","lam_interceptsb","sigma_interceptsb",
                        "sigma_y"),"mean"]; t(param)

