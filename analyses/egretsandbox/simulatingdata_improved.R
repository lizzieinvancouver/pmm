# Feb 2, 2025
# Started by V Van der Meersch
# Working on simulating data - closer to real data we couled get?
# (now in function simulate_data_improved, in toolbox)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("lizzie", getwd())>0)) {
  setwd("~/Documents/git/teaching/stan/pmm/analyses")
} else if(length(grep("victor", getwd())>0)) {
  setwd("/home/victor/projects/pmm/analyses/egretsandbox")
} else{
  setwd("/home/deirdre/pmm") # for midge
}

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)
require(shinystan)
require(future.apply)
require(progressr)

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
handlers(global = TRUE)
handlers("txtprogressbar")
set.seed(1234567)

# parameters
params <- list(a_z = rnorm(1,0,1.5), # root value intercept
               lambda_a = rbeta(1,1.5,1.5), # lambda intercept
               sigma_a = abs(rnorm(1,0,1)), # rate of evolution intercept
               b_z = rnorm(1,0.5,1), # root value trait1 slope
               lambda_b = rbeta(1,1.5,1.5), # lambda trait1
               sigma_b = abs(rnorm(1,0,1)) # rate of evolution trait1
)

## Simulate some data ## 
# phylogenetic structure
nspecies <- 20 # number of species
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1) # pure birth model
spetree$tip.label <- paste("s", 1:nspecies, sep="")
scaledtree_intercept <- rescale(spetree, model = "lambda", params$lambda_a)         
params$intercepts <- fastBM(scaledtree_intercept, a = params$a_z, mu = 0, sig2 = params$sigma_a ^ 2)
scaledtree_slope <- rescale(spetree, model = "lambda", params$lambda_b)         
params$slopes <- fastBM(scaledtree_slope, a = params$b_z, mu = 0, sig2 = params$sigma_b ^ 2)

# we consider that a study looks at only one species, but a species can be considered by several studies
nstudies_perspecies <- round(runif(nspecies, 1, 4)) # number of different studies per species
nstudies <- sum(nstudies_perspecies) # total number of studies
species_study <- rep(1:nspecies, times = nstudies_perspecies) # species id considered by each study
ntreat_perstudies <- round(runif(nstudies, 2, 6)) # number of different treatments applied per study
nexps <- sum(ntreat_perstudies) # total number of experiments (i.e. studies*treatments)
ntrialseeds_perexp <- round(runif(nexps, 10, 100)) # number of seeds per experiments
species <- rep(species_study, times = ntreat_perstudies)

# experimental observations
x <- rnorm(n = nexps, mean = 0, sd = 1) # treatments applied
yhat <- plogis(intercepts[species] + x * slopes[species], location = 0, scale = 1)
y = rbinom(n = nexps, size = ntrialseeds_perexp, prob = yhat)

plot(y/ntrialseeds_perexp ~ x, col = species)
