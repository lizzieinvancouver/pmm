## Started 2 August 2022 ##
## By Lizzie and Dinara Mamatova ##
## Copy of ubermini_2.R for testing
## Generate test data for phylogeny Stan models ##

## ADD intercept to estimate
# y ~ a[phylo] +  b_force[phylo]*x + error

set.seed(2021)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("dindin", getwd())>0)) { 
  setwd("/Users/dindin/Documents/work/pmm/analyses") 
} else if(length(grep("mamatova", getwd())>0)) {
  setwd("/home/mamatova/speedtests") 
} else setwd("~/Documents/git/teaching/stan/pmm/analyses")


## load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

#options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

nspecies = 50
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait
m <- 0.6  # root trait value, compare to b_z (set to -0.4 if you want to recreate early Nov 2020 issue in OSPREE data: `multi_normal_lpdf: LDLT_Factor of covariance parameter is not positive definite'); as of late July 2022; changing this to -0.6 returns good values
lam <- 0.7 # lambda
sigma_y <- 0.01 # sigma_y
sigma_interceptsb <- 0.1 # rate of evolution (a constant rate across the tree)

scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sigma_interceptsb ^ 2)
phylosig(x=slopez, tree=spetree, method="lambda")

# ADAPT here to include intercept also...
# set up the intercept
a_z = 4 # root value intercept
lam_interceptsa = 0.4 # lambda intercept
sigma_interceptsa = 0.2 # rate of evolution intercept

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", lam_interceptsa)         
intercepts <- fastBM(scaledtree_intercept, a = a_z, mu = 0, sig2 = sigma_interceptsa ^ 2)
phylosig(x=intercepts, tree=scaledtree_intercept, method="lambda")

# Set up priors
priors <- list(
  
  # a_z_prior_mu = 4, # true value # true val of intercept 
  # a_z_prior_sigma = 1, # root value intercept sigma
  # lam_interceptsa_prior_alpha = 4, # lambda intercept prior alpha
  # lam_interceptsa_prior_beta = 6, # lambda intercept prior beta
  # sigma_interceptsa_prior_mu = 0.2, # true value of mu for the rate of evolution intercept
  # sigma_interceptsa_prior_sigma = 0.2, # rate of evloution prior 
  # b_z_prior_mu = 0.6, # true value
  # b_z_prior_sigma = 1,
  # lam_interceptsb_prior_alpha = 7, #
  # lam_interceptsb_prior_beta = 3, # 
  # sigma_interceptsb_prior_mu = 0.1, # true value
  # sigma_interceptsb_prior_sigma = 0.1,
  # sigma_y_mu_prior = 0.01,
  # sigma_y_mu_sigma = 1
  
  a_z_prior_mu = 4,
  a_z_prior_sigma = 5,
  lam_interceptsa_prior_alpha = 2,
  lam_interceptsa_prior_beta = 2,
  sigma_interceptsa_prior_mu = 0.2,
  sigma_interceptsa_prior_sigma = 1,
  b_z_prior_mu = 0.6,
  b_z_prior_sigma = 1,
  lam_interceptsb_prior_alpha = 2,
  lam_interceptsb_prior_beta = 2,
  sigma_interceptsb_prior_mu = 0.1,
  sigma_interceptsb_prior_sigma = 1,
  sigma_y_mu_prior = 0.01,
  sigma_y_mu_sigma = 0.5
)

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

dfhere <- data.frame(x=numeric(), y=numeric())

for (i in 1:length(slopez)){
    slopehere <- slopez[i]
    intercpthere <- intercepts[i]
    xhere <- rnorm(nind, 10, 3) # these are experiments, no phylo structure in x 
    yhere <- intercpthere + xhere*slopehere # ADAPT here to include intercept
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

start_time <- Sys.time()
# Fit model
testme <- stan("stan/uberspeed.stan", # Or ubermini_2.stan
               data = append(list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                          x=dfhere$x, y=dfhere$yerr,
                         Vphy=vcv(spetree, corr = TRUE)), # Note: In this case corr=TRUE/FALSE produce same output
                         priors), 
               iter=4000, chains=4,
               seed = 123456,
# Only used to fix the chain divergences
               control = list(adapt_delta = 0.95)) # seed=123456
end_time <- Sys.time()

# Find time difference
end_time - start_time

# Summarize fit
summary(testme)$summary
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"mean"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"50%"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"25%"]
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z", "sigma_y"),"75%"]

#pairs(testme, pars = c("a_z","lam_interceptsa","sigma_interceptsa", "b_z","lam_interceptsb","sigma_interceptsb","sigma_y", "lp__")) 

fitContinuous(spetree, slopez, model="lambda")
fitContinuous(spetree, intercepts, model="lambda")

# Compare true slopes to estimated slopes
sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)), "mean"]
# slopez

plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])
abline(a = 0, b = 1, lty = "dotted")

plot(intercepts, sumer[grep("a\\[(.*?)\\]", rownames(sumer)), "mean"])
abline(a = 0, b = 1, lty = "dotted")

