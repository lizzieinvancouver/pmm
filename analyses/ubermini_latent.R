## Started mid October 2020 ##
## By Geoff ##
## Generate test data for my phylogeny Stan models ##

## Someday we want:
# y ~ a[sp] + b_force[phylo]*x1 + b_chill[phylo]*x2 + b_photo[phylo]*x3 + error

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")

## Load packages
library(ape)
library(geiger)
library(phytools)
library(rstan)
library(shinystan)

## Set number of cores to use
options(mc.cores = 4)
## options(browser = "chromium")

## set seed
set.seed(20201008)
## Parameters
nspecies = 80
nind = 15
### Trait stuff
param <- list(
    b_force_latent_mu = 2.7,
    b_force_latent_sigma = .1,
    sigma_y = .8,
    sigma_phylo = 0.3)
## Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")
### Generate phylogenetic effect using (unscaled) tree
slope.phylogeny <- fastBM(spetree, a = 0, mu = 0, sig2 = param[["sigma_phylo"]] ^ 2) # I am not sure if sig2 is really what we want to vary here, but maybe you have read up on it and it is ... a = 0 sets the root trait value to 0, this might be what we want given your 'actual slope' additive model? 
### Generate latent slope
slope.latent <- rnorm(n = nspecies, mean = param[["b_force_latent_mu"]], sd = param[["b_force_latent_sigma"]])
### Add together
slope.actual <- slope.phylogeny + slope.latent 

## What is lambda? Model will struggle if too low (because it will want sigma_phylo to be 0, the lower bound)
phylosig(x=slope.actual, tree=spetree, method="lambda")$lambda

## Organize into data frame
dfhere <- data.frame(x=numeric(), y=numeric())
for (i in 1:length(slope.actual)){
    slopehere <- slope.actual[i]
    xhere <- rnorm(nind, 10, 3)
    yhere <- xhere*slopehere
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}
dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, param[["sigma_y"]])

## Fit model
fit1 <- stan("stan/pmm_geoff.stan",
             data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                       x=dfhere$x, y=dfhere$yerr,
                       Vphy=vcv(spetree, corr = TRUE)),
             iter=5000, warmup = 4000, chains=4, seed=123456)

## Compare true values to estimates
### True values
summary(fit1, pars = list("b_force_latent_mu", "b_force_latent_sigma", "sigma_y", "sigma_phylogeny"))$summary
### Estimates
param

## Compare lambda estimates
### Based on tree
phylosig(x=slope.actual, tree=spetree, method="lambda")$lambda
### Based on fits
fit1.extract <- extract(fit1)
(mean(fit1.extract$sigma_phylogeny) ^ 2) / ((mean(fit1.extract$b_force_latent_sigma) ^ 2) + (mean(fit1.extract$sigma_phylogeny) ^ 2))

## View diagnostics
launch_shinystan(fit1)
