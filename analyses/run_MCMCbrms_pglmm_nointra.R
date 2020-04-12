## Code from Joly et al. 2019 ##
## He sent this in an email on 10 March 2020 ##
## Edits by Nacho & Lizzie so far ##

## This is very similar to run_MCMCbrms_pglmm.R but the simulation removes intraspecific structure ##
## I cannot get the PGLS models to run so have removed all by the GLS from here ##
## See pmm-pgls-simulations_v2.R for some of the email exchange with Joly about this ... ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")


## load packages
library(shinystan)
library(caper)
library(brms)
library(pez)
library(rstan)
library(phytools)
library(MCMCglmm)
library(dplyr)
library(knitr)
library(broom)
library(loo)
# Packages required
require(nlme)
require(Matrix)
require(phyclust)


## source aux functions
source("source/corIntra.R")
source("source/pmm-get-simulated-data-nointra.R")


## Parameters (set your own set of parameters according to PMM_models_comparison.csv)
# nspecies: number of species in the tree
# nindividuals: number of individuals sampled per species
# sigma.sq.x: Rate of evolution in x
# sigma.sq.p: Phylogenetic variance
# sigma.sq.c: Intraspecific correlation structure variance
# sigma.sq.e: Rate of residual variation
# B: Regression slope
nspecies=20
nindividuals=20
B = 0.75
sigma.sq.x = 2
sigma.sq.p = 1
sigma.sq.c = 1
sigma.sq.e = 1


## get data using your specified parameters
# note: simdata (aka, now simd below) returns only simdata$data on FIRST run, run it again and you get the full list specificed in f(x) ... or you get the full list every time if you don't call it simdata (something to do with <<- in the f(x) code?!)

simd <- one.sim.pmm.nointra(nspecies = nspecies, nindividuals = nindividuals,B = B,
                     sigma.sq.x = sigma.sq.x,sigma.sq.p = sigma.sq.p,
                     sigma.sq.e = sigma.sq.e)

  
#plot(simd$data$x,simd$data$y)
#simd <- simd.nointra

## run MCMCglmm models 
    
###
# Run the MCMCpglmm analyses with different correlation structures.
# For each, we will use a very diffuse prior.
    
# Model M.0 -> no random effects
    priorpr.m0 <- list(R = list(V = 1, nu = 0.002))
    M.0 <- MCMCglmm(y ~ x, 
                    data=simd$data, scale=TRUE,prior=priorpr.m0,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
    
# Model M.1 -> adding phylogenetic random effects
    priorpr.m1 <- list(R = list(V = 1, nu = 0.002), 
                       G = list(G1 = list(V = 1, nu = 0.002)))
    M.1 <- MCMCglmm(y ~ x, random = ~ animal, pedigree = simd$phylo,
                    data=simd$data, scale=TRUE,prior=priorpr.m1,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)

## Fit GLS
#
gls.mod0 <- gls.mod1 <- gls.mod2 <- gls.mod3 <-NULL
try(gls.mod0 <- gls(y~x,data=simd$data)) # no phylogenetic correlation

## Add brms (first we do the weird within-group centering, wgc)
# 
simd$data$specmean <- 
  with(simd$data, sapply(split(simd$data$x, simd$data$animal), mean)[simd$data$animal])
simd$data$specphylo <- simd$data$animal
simd$data$withspecmean <- simd$data$x-simd$data$specmean

A <- ape::vcv.phylo(simd$phylo)

brmmod <- brm(
  y ~ x + (1|animal) + (1|specphylo), 
  data = simd$data, family = gaussian(), 
  cov_ranef = list(phylo = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 4, cores = 4, 
  iter = 4000, warmup = 1000
)


brmmod.wgc <- brm(
  y ~ specmean + withspecmean + (1|animal) + (1|specphylo), 
  data = simd$data, family = gaussian(), 
  cov_ranef = list(phylo = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 4, cores = 4, 
  iter = 4000, warmup = 1000
)


hyp <- paste(
  "sd_specphylo__Intercept^2 /", 
  "(sd_specphylo__Intercept^2 + sd_animal__Intercept^2 + sigma^2) = 0")
(hyp <- hypothesis(brmmod, hyp, class = NULL))

hypwgc <- paste(
  "sd_specphylo__Intercept^2 /", 
  "(sd_specphylo__Intercept^2 + sd_animal__Intercept^2 + sigma^2) = 0")
(hypwgc <- hypothesis(brmmod.wgc, hypwgc, class = NULL))


###
# Compile the results of PMM (in MCMCglmm) and PGLS
  
# Model summaries
  M.0.sum <- summary(M.0)
  M.1.sum <- summary(M.1)
  gls.mod0.sum <- summary(gls.mod0)
  gls.mod1.sum <- summary(gls.mod1)
  brmmod.sum <- summary(brmmod)
  brmmod.wgc.sum <- summary(brmmod.wgc)
  
# Data frame
# 
  
  res <- data.frame(models=c("M.0","M.1", "gls.mod0", "brms", "brms.wgc"),
                    #random.effects=c("NA","inter","infra","inter+infra","ols","glsInter","glsIntra","glsInterIntra"),
                    slope=c(M.0.sum$solutions[2,1],M.1.sum$solutions[2,1],
                            gls.mod0$coefficients[2],
                            brmmod.sum$fixed[2,1], brmmod.wgc.sum$fixed[3,1]),
                    slope.sdse=c(sd(M.0$Sol[,'x']),sd(M.1$Sol[,'x']),
                               gls.mod0.sum$tTable[2,2],
                               brmmod.sum$fixed[2,2],
                               brmmod.wgc.sum$fixed[3,2]))
  
  # Calculate the absolute difference with the 'true' slope
  res$slope.diff <- abs(B-res$slope)
  res$best.slope <- ifelse(res$slope.diff==min(res$slope.diff),1,0)
  
  # Calculate the variance contributions
  var.contr <- function(model){
    xx <- model$VCV/apply(model$VCV,1,sum)
    apply(xx,2,mean)
  }
  M.0.var <- var.contr(M.0)
  M.1.var <- var.contr(M.1)

  # Extra phylo stuff
  res$Rstruct <- c(M.0.sum$Rcovariances[1,1],M.1.sum$Rcovariances[1,1],
               NA, NA, NA)
  res$Gstruct <- c(NA,M.1.sum$Gcovariances[1,1],
               NA, NA, NA)
  res$brmsphylo <- c(NA,NA, NA, 
              brmmod.sum$random$specphylo[1,1], brmmod.wgc.sum$random$specphylo[1,1])
  res$brmsanimal <- c(NA,NA, NA, 
             brmmod.sum$random$animal[1,1], brmmod.wgc.sum$random$animal[1,1])

  # Heritability (h^2)
  res$h2 <- c(NA,M.1.var["animal"]/sum(M.1.var), NA,
              hyp$hypothesis$Estimate, hypwgc$hypothesis$Estimate)

  # Lizzie thinks h2 should equal (which seems vaguely correct based on a few tests using MCMCglmm)
  sigma.sq.p/(sigma.sq.p+sigma.sq.e) 
  
# write.csv(res,file = "output/sim1_IMC.csv")


  
# END
