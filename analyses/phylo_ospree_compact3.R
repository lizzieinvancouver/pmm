## Started 28 November 2020 ##
## By Lizzie ## 
## Taken mainly from Phylo_ospree_reanalyses.R ##

library("ape")
library("rstan")

options(mc.cores = 4)

# Data
phylo <- read.tree("input/ospreephylo.tre")
bb.stan <- read.csv("input/ospreebb.csv")

# Step 1: Get spps and VCVPHY in same order
# bb.stan$spps[phylo$tip.label]
phylo$tip.label
d <- bb.stan[match(phylo$tip.label, bb.stan$spps),]

phymatch <- data.frame(tip=phylo$tip.label, sppnum=c(1:length(phylo$tip.label)))
d <- merge(bb.stan, phymatch, by.x="spps", by.y="tip")
d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)

## Two slope, one intercept model
### Set priors
phypriors <- list(
    a_z_prior_mu = 30,
    a_z_prior_sigma = 5,
    b_zf_prior_mu = -4,
    b_zf_prior_sigma = 5,
    b_zc_prior_mu = -8,
    b_zc_prior_sigma = 5,
    b_zp_prior_mu = -3,
    b_zp_prior_sigma = 5,
    sigma_interceptsa_prior_mu = 40,
    sigma_interceptsa_prior_sigma = 5,
    sigma_interceptsbf_prior_mu = 5,
    sigma_interceptsbf_prior_sigma = 5,
    sigma_interceptsbc_prior_mu = 5,
    sigma_interceptsbc_prior_sigma = 5,
    sigma_interceptsbp_prior_mu = 5,
    sigma_interceptsbp_prior_sigma = 5,
    sigma_y_mu_prior = 20,
    sigma_y_mu_sigma = 1)

### Fit model
testme <- stan("stan/uber_threeslopeintercept.stan",
    data=append(list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x1=d$force.z, x2 = d$chill.z, x3=d$photo.z,
       y=d$resp, Vphy=vcv(phylo, corr = TRUE)), phypriors),
    iter=2000, chains=4,
    refresh = 10)

saveRDS(testme, "compact3.rds")

## ### Add species names
## names(testme)[4:(nspecies+3)] <- unique(d$phylo)

## ### Summarize full fit
## summary(testme)$summary
        
## ### Summarize lambda, b_zf, b_zc, intercepts, and sigmas
## summary(testme, pars = list("lam_interceptsa", "lam_interceptsb", "sigma_interceptsa", "sigma_interceptsbf", "sigma_interceptbc", "a_z", "b_zf", "b_zc" "sigma_y"))$summary


