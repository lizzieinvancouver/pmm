## Started 29 May 2020 ##
## Generate test data for my phylogeny Stan models ##

## Someday we want:
# y ~ a[sp] + b_force[phylo]*x1 + b_chill[phylo]*x2 + b_photo[phylo]*x3 + error

## Uber simple to start ...
## skip an intercept and just estimate
# y ~ b_force[phylo]*x1 + error

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")

## load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)
library(truncnorm)

options(mc.cores = parallel::detectCores())
## options(mc.cores = 4)

nspecies = 90
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait
m <- 0.6  # root trait value, compare to b_z (set to -0.4 if you want to recreate early Nov 2020 issue in OSPREE data: `multi_normal_lpdf: LDLT_Factor of covariance parameter is not positive definite')
lam <- 0.7 # lambda
sigma_y <- 0.01 # sigma_y
sigma_interceptsb <- 0.1 # rate of evolution (a constant rate across the tree)

scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sigma_interceptsb ^ 2)
phylosig(x=slopez, tree=spetree, method="lambda")

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

dfhere <- data.frame(x=numeric(), y=numeric())

for (i in 1:length(slopez)){
    slopehere <- slopez[i]
    xhere <- rnorm(nind, 10, 3) # these are experiments, no phylo structure in x 
    yhere <- xhere*slopehere
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}

dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigma_y)

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    b_z.temp <- rnorm(n = 1, mean = 0, sd = 3)
    sigma_interceptsb.temp <- runif(n = 1, min = 0, max = 1)
    return(list("sigma_y" = runif(n = 1, min = 0, max = 1),
                "lam_interceptsb" = rbeta(n = 1, shape1 = 2, shape2 = 2),
                "sigma_interceptsb" = sigma_interceptsb.temp,
                "b_z" = b_z.temp,
                "b_force" = rnorm(n = nspecies, mean = b_z.temp, sd = sigma_interceptsb.temp)))
}

# Fit model
testme <- stan("stan/ubermini_2_biggerpriors.stan", # Or ubermini_2.stan
               data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                          x=dfhere$x, y=dfhere$yerr,
                         Vphy=vcv(spetree, corr = TRUE)), # Note: In this case corr=TRUE/FALSE produce same output
               iter=2000, chains=4,
               init = simu_inits
               ) # seed=123456
# Summarize fit
summary(testme)$summary
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"mean"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"50%"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"25%"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"75%"]

fitContinuous(spetree, slopez, model="lambda")

# Compare true slopes to estimated slopes
sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)), "mean"]
# slopez

plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])
abline(a = 0, b = 1, lty = "dotted")
