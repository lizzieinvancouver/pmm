## Started 28 November 2020 ##
## By Lizzie ## 
## Taken mainly from Phylo_ospree_reanalyses.R ##


## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory.
# Add in your own path in an if statement for your file structure (please don't delete this)
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/teaching/stan/pmm/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

# libraries
library("ape")
library("rstan")

options(mc.cores = parallel::detectCores())
## options(mc.cores = 4)

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
    a_z_prior_sigma = 10,
    lam_interceptsa_prior_alpha = 1, # 
    lam_interceptsa_prior_beta = 1, # 
    sigma_interceptsa_prior_mu = 40,
    sigma_interceptsa_prior_sigma = 20,
    b_zf_prior_mu = -4,
    b_zf_prior_sigma = 5,
    lam_interceptsbf_prior_alpha = 1, #
    lam_interceptsbf_prior_beta = 1, # 
    sigma_interceptsbf_prior_mu = 5,
    sigma_interceptsbf_prior_sigma = 5,
    b_zc_prior_mu = -8,
    b_zc_prior_sigma = 5,
    lam_interceptsbc_prior_alpha = 1, #
    lam_interceptsbc_prior_beta = 1, # 
    sigma_interceptsbc_prior_mu = 5,
    sigma_interceptsbc_prior_sigma = 5,
    b_zp_prior_mu = -3,
    b_zp_prior_sigma = 5,
    lam_interceptsbp_prior_alpha = 1, #
    lam_interceptsbp_prior_beta = 1, # 
    sigma_interceptsbp_prior_mu = 5,
    sigma_interceptsbp_prior_sigma = 5,
    sigma_y_mu_prior = 20,
    sigma_y_mu_sigma = 1)

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
    b_zf.temp <- rnorm(n = nspecies, mean = phypriors[["b_zf_prior_mu"]], sd = phypriors[["b_zf_prior_sigma"]])
    b_zc.temp <- rnorm(n = nspecies, mean = phypriors[["b_zc_prior_mu"]], sd = phypriors[["b_zc_prior_sigma"]])
    b_zp.temp <- rnorm(n = nspecies, mean = phypriors[["b_zp_prior_mu"]], sd = phypriors[["b_zp_prior_sigma"]])
    return(append(list(a = a_z.temp,
                       b_force = b_zf.temp,
                       b_chill = b_zc.temp,
                       b_photo = b_zp.temp),
                  phypriors))
}

## Fit model
testme <- stan("stan/uber_threeslopeintercept.stan",
               data=append(list(N=nrow(d),
                                n_sp=nspecies,
                                sp=d$sppnum,
                                x1=d$force.z,
                                x2 = d$chill.z,
                                x3=d$photo.z,
                                y=d$resp,
                                Vphy=vcv(phylo, corr = TRUE)),
                           phypriors),
               init = simu_inits,
               iter = 6000,
               warmup = 5000,
               chains = 4,
               seed = 62921,
               refresh = 200
               )

## Save fitted posterior
saveRDS(testme, "output/testme.rds")

## Summarize full fit
summary(testme)$summary

## Summarize lambdas, b_zf, b_zc, , b_zp, intercept mean, and sigmas
summary(testme, pars = list("a_z", "lam_interceptsa", "sigma_interceptsa", "b_zf", "lam_interceptsbf", "sigma_interceptsbf", "b_zc", "lam_interceptsbc", "sigma_interceptsbc", "b_zp", "lam_interceptsbp", "sigma_interceptsbp", "sigma_y"))$summary

# Reinstating some useful plotting code
# See also https://github.com/lizzieinvancouver/pmm/blob/5014539f8a7cfc659298d20d49a0935a8ced305d/analyses/phlyo_opsree_compact.R
source("source/stan_utility.R")
fit <- readRDS("output/testme.rds")
check_all_diagnostics(fit)

names(fit)[grep(pattern = "^a\\[", x = names(fit))] <- phylo$tip.label
names(fit)[grep(pattern = "^b_force", x = names(fit))] <- phylo$tip.label
names(fit)[grep(pattern = "^b_chill", x = names(fit))] <- phylo$tip.label
names(fit)[grep(pattern = "^b_photo", x = names(fit))] <- phylo$tip.label

pdf(file = "output/estimates1.pdf", onefile = TRUE, height = 35, width = 6)
plot(fit, pars = c("a_z", "a"))
plot(fit, pars = c("b_zf", "b_force"))
plot(fit, pars = c("b_zc", "b_chill"))
plot(fit, pars = c("b_zp", "b_photo"))
dev.off()
pdf(file = "output/estimates2.pdf", onefile = TRUE, height = 11.5, width = 8)
plot(fit, pars = c("lam_interceptsa", "lam_interceptsbf", "lam_interceptsbc", "lam_interceptsbp"))
plot(fit, pars = c("sigma_interceptsa", "sigma_interceptsbf", "sigma_interceptsbc", "sigma_interceptsbp", "sigma_y"))
dev.off()





fitsamples <- extract(fit)
partition <- partition_div(fit)
div_samples <- partition[[1]]
nondiv_samples <- partition[[2]]

c_dark <- c("#8F272780")
green <- c("#00FF0080")

partition <- partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

pars <- c("sigma_y", "lam_interceptsa", "sigma_interceptsa", "lam_interceptsbf", "sigma_interceptsbf",
    "lam_interceptsbc", "sigma_interceptsbc", "lam_interceptsbp", "sigma_interceptsbp", "b_zf", "b_zc",
    "b_zp", "a_z", "lp__")
sumfit <- summary(fit)$summary
sumfit[which(rownames(sumfit) %in% pars),]

colz <- c(rep("#8F272780", nrow(nondiv_params[pars])), rep("#00FF0080", nrow(div_params[pars])))

pdf(file="figures/ospreepairs.pdf")
pairs(~sigma_y + lam_interceptsa + sigma_interceptsa + lam_interceptsbf +  sigma_interceptsbf + 
    lam_interceptsbc +  sigma_interceptsbc +  lam_interceptsbp +  sigma_interceptsbp + 
    b_zf + b_zc + b_zp + a_z + lp__, 
    data=rbind(nondiv_params[pars], div_params[pars]), pch=21, cex=0.8, bg=colz)
dev.off()


