## Started 15 May 2022 ##
## Adapted from simulate_fit_oneslopeintercept_cholesky_dl
## By Lizzie ##
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

# Set seed
# set.seed(2021)

nfam = 50
nspeciesinfam = 10
nind = 10
nspecies = nfam*nspeciesinfam


# Simulate species tree with a pure birth model
spetree <- pbtree(n=nfam, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nfam, sep="")

# Now set up the trait parameters
param <- list(a_z = 4, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              sigma_y = 0.01 # overall sigma
              )

speciesinfamsig = 1.5

# Set priors
phypriors <- list(
    a_z_prior_mu = 4, # true value
    a_z_prior_sigma = 3,
    lam_interceptsa_prior_alpha = 4, # 
    lam_interceptsa_prior_beta = 6, # 
    sigma_interceptsa_prior_mu = 0.2, # true value
    sigma_interceptsa_prior_sigma = 0.2,
    sigma_y_prior_mu = 0.01,
    sigma_y_prior_sigma = 1  
)

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)



# Generate tree, same as before, but we'll call it the family-level tree
dfhere <- data.frame(fam = c(), sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
for (i in 1:nfam){
    famhere <- rnorm(nspeciesinfam, 0, speciesinfamsig)
    temp <- data.frame(fam = rep(i, nspeciesinfam*nind),
                       sp = rep(1:nspeciesinfam, each=nind),
                       intercept = rep(intercepts[i], nind),
                       spadd = rep(famhere, each=nind))
    dfhere <- rbind(dfhere, temp)
}

# set the intercept to zero to make my life easier ... 
dfhere$mu <- dfhere$intercept + dfhere$spadd
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

# make species unique (not always critical to do, depending on model structure, but safer)
unique(dfhere$sp)
dfhere$famsp <- as.numeric(as.factor(paste(dfhere$fam, dfhere$sp, sep="")))
unique(dfhere$famsp)


# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = 1)
    return(append(list(a = a_z.temp), param))
}

# Fit model

test <- stan("sandbox/stan/oneinterceptcholforsyncfam.stan",
               data = append(list(N = nrow(dfhere),
                                  Nfam = nfam,
                                  family = dfhere$fam,
                                  Nspp = nspecies,
                                  species = dfhere$famsp,
                                  ypred = dfhere$y,
                                  Vphy = vcv(spetree, corr = TRUE)),
                             phypriors),
               # init = simu_inits,
               iter = 3000,
               warmup = 2000,
               chains = 4,
               seed = 62921
             )


summary(test)$summary[c("a_z","lam_interceptsa","sigma_interceptsa",
                        "sigma_y"),"mean"]; t(param)

# I would like to figure out how to test the lower-level parameters ...
# Can't tell if this is awful, or it's the cholesky
sumer <- summary(test)$summary
plot(intercepts-4, sumer[grep("a\\[", rownames(sumer)), "mean"])
abline(a = 0, b = 1, lty = "dotted")
# Really small values...
hist(sumer[grep("a\\[", rownames(sumer)), "mean"])
hist(sumer[grep("tau_family\\[", rownames(sumer)), "mean"])


