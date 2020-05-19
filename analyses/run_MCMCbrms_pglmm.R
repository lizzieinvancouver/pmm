## Code from Joly et al. 2019 ##
## See ee pmm-pgls-simulations_v2.R for email with Joly ##
## He sent this in an email on 10 March 2020 ##
## Edits by Nacho & Lizzie so far ##

## Modified by Nacho and Lizzie 8-12 April - to source functions to simulate data and run:
# MCMCglmm
# PGLS
# BRMS models
## 

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else if
(length(grep("catchamberlain", getwd()))>0) {setwd("~/Documents/git/pmm/analyses")
}  else
setwd("~/Documents/git/teaching/stan/pmm/analyses")

set.seed(113)

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
source("source/pmm-get-simulated-data.R")
# source("source/pmm-get-simulated-data-nointra.R")
# source("source/pmm-get-simulated-data-joly.R")


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
sigma.sq.x = 1
sigma.sq.p = 0.25
sigma.sq.c = 1
sigma.sq.e = 1
scaletree = 1
phyloSlope = FALSE


## get data using your specified parameters
# note: simdata (aka, now simd below) returns only simdata$data on FIRST run, run it again and you get the full list specificed in f(x) ... or you get the full list every time if you don't call it simdata (something to do with <<- in the f(x) code?!)

simd <-one.sim.pmm(nspecies = nspecies, nindividuals = nindividuals,B = B,
                     scaletree=scaletree,
                     sigma.sq.x = sigma.sq.x,sigma.sq.p = sigma.sq.p,
                     sigma.sq.c = sigma.sq.c,sigma.sq.e = sigma.sq.e,
                     phyloSlope=phyloSlope)

  
# plot(simd$data$x,simd$data$y)

## run MCMCglmm models 
# Prepare infraspecific correlation matrix for MCMCglmm
#
# Single value decomposition of the intraspecific structure
#   -> follow code of Stone et al. 2012
    intra.svd <- svd(simd$intramat)
    intra.svd <- intra.svd$v %*% (t(intra.svd$u) * sqrt(intra.svd$d))
    rownames(intra.svd) <- colnames(intra.svd) <- rownames(simd$intramat)

# The following is important and MCMCglmm searches the variables in the 
# global environment
    intra.svd <<- intra.svd
    
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


if(FALSE){ # Lizzie's poor attempts at writing my own stan code, need to check Rethinking next, I think
# get numeric species names in the same order as the tip labels... 
simd$data$stansp <- paste("s", as.numeric(as.factor(factor(simd$data$animal, levels = unique(simd$data$animal)))), sep="")
simd$data$stansp_num <- as.numeric(gsub('s', '', simd$data$stansp))
    
testme <- stan("stan/nointer_2level_forceTEST.stan",
                data=list(N=nrow(simd$data), n_sp=nspecies, sp=simd$data$stansp_num,
                force=simd$data$x, y=simd$data$y,
                Vphy=vcv(simd$phylo)),
                iter=1000, chains=4, seed=123456)
summary(testme)$summary

# From Will's code
# ...the estimated phylogenetic signal (Pagel's Lambda) of the
# estimated environmental response is ~.7 (i.e., there's signal)
# (see Hadfield & Nakagawa (2010) J. Evol. Biol. 23(3): 494-508 for
# this method of calculating it)
lam.int <- mean(extract(testme)[["lam_intercepts"]])
null.int <- mean(extract(testme)[["null_intercepts"]])
lam.int / (null.int + lam.int) # very wrong ... I believe (says Lizzie)
    }
    
## run PGLS models 
#

inter.mat.inflated <- inflate.mat(simd$intermat,nspecies,nindividuals)

# Generate phylogenetic correlation matrices
# The intra correlation structure takes the inter-specific and the intra-
# specific correlation matrices as input.
pgls.intra.cor <- corIntra(0, inter=inter.mat.inflated, intra=simd$intramat, fixed=TRUE)
pgls.inter.cor <- corIntra(1, inter=inter.mat.inflated, intra=simd$intramat, fixed=TRUE)
pgls.intrainter.cor <- corIntra(0.5, inter=inter.mat.inflated, intra=simd$intramat, fixed=FALSE)
    
# Fit PGLS
gls.mod0 <- gls.mod1 <- gls.mod2 <- gls.mod3 <-NULL
try(gls.mod0 <- gls(y~x,data=simd$data)) # no phylogenetic correlation
try(gls.mod1 <- gls(y~x,correlation = pgls.inter.cor, data=simd$data),silent=TRUE) # inter structure 
try(gls.mod2 <- gls(y~x,correlation = pgls.intra.cor, data=simd$data),silent=TRUE) # intra structure
try(gls.mod3 <- gls(y~x,correlation = pgls.intrainter.cor, data=simd$data),silent=TRUE) # Estimating delta


## Add brms (last we do the weird within-group centering, wgc)
# 
simd$data$specmean <- 
  with(simd$data, sapply(split(simd$data$x, simd$data$animal), mean)[simd$data$animal])
simd$data$withspecmean <- simd$data$x-simd$data$specmean
simd$data$specphylo <- simd$data$animal # this name needs to match up in cov_ranef() below

A <- ape::vcv.phylo(simd$phylo)

brmmod <- brm(
  y ~ x + (1|specphylo), 
  data = simd$data, family = gaussian(), 
  cov_ranef = list(specphylo = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 4, cores = 4, 
  iter = 4000, warmup = 1000
)


brmmod.simple <- brm(
  y ~ x + (1|animal) + (1|specphylo), 
  data = simd$data, family = gaussian(), 
  cov_ranef = list(specphylo = A),
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
  cov_ranef = list(specphylo = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 4, cores = 4, 
  iter = 4000, warmup = 1000
)

# See also: https://discourse.mc-stan.org/t/phylogenetic-model-with-repeated-measurements-without-using-species-mean/4801

hyp <- paste(
  "sd_specphylo__Intercept^2 /", 
  "(sd_specphylo__Intercept^2 + sigma^2) = 0")
(hyp <- hypothesis(brmmod, hyp, class = NULL))


hypsimple <- paste(
  "sd_specphylo__Intercept^2 /", 
  "(sd_specphylo__Intercept^2 + sd_animal__Intercept^2 + sigma^2) = 0")
(hypsimple <- hypothesis(brmmod.simple, hypsimple, class = NULL))


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
  gls.mod2.sum <- summary(gls.mod2)
  gls.mod3.sum <- summary(gls.mod3)
  brmmod.sum <- summary(brmmod)
  brmmod.simple.sum <- summary(brmmod.simple)
  brmmod.wgc.sum <- summary(brmmod.wgc)
  
# Data frame
# 
  
  res <- data.frame(models=c("M.0","M.1","gls0","gls1","gls2","gls3", "brms.1ranef", "brms.2ranef", "brms.wgc"),
                    #random.effects=c("NA","inter","infra","inter+infra","ols","glsInter","glsIntra","glsInterIntra"),
                    slope=c(M.0.sum$solutions[2,1],M.1.sum$solutions[2,1],
                            gls.mod0$coefficients[2],gls.mod1$coefficients[2],
                            gls.mod2$coefficients[2],gls.mod3$coefficients[2],
                            brmmod.sum$fixed[2,1],  brmmod.simple.sum$fixed[2,1],
                            brmmod.wgc.sum$fixed[3,1]),
                    slope.sdse=c(sd(M.0$Sol[,'x']),sd(M.1$Sol[,'x']),
                               gls.mod0.sum$tTable[2,2],
                               gls.mod1.sum$tTable[2,2],
                               gls.mod2.sum$tTable[2,2],
                               gls.mod3.sum$tTable[2,2],
                               brmmod.sum$fixed[2,2],
                               brmmod.simple.sum$fixed[2,2],
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
              NA,NA,NA,NA, NA, NA, NA)
  res$Gstruct <- c(NA,M.1.sum$Gcovariances[1,1],
              NA,NA,NA,NA, NA, NA, NA)
  res$brmsphylo <- c(NA,NA,
              NA,NA,NA,NA,brmmod.sum$random$specphylo[1,1],
              brmmod.simple.sum$random$specphylo[1,1], brmmod.wgc.sum$random$specphylo[1,1])
  res$brmsanimal <- c(NA,NA,
              NA,NA,NA,NA, NA, brmmod.simple.sum$random$animal[1,1], brmmod.wgc.sum$random$animal[1,1])

  # Heritability (h^2)
  res$h2 <- c(NA,M.1.var["animal"]/sum(M.1.var),
              NA,NA,NA,NA, hyp$hypothesis$Estimate,
              hypsimple$hypothesis$Estimate, 
              hypwgc$hypothesis$Estimate)

  # Lizzie thinks h2 should equal:
  sigma.sq.p/(sigma.sq.p+sigma.sq.c+sigma.sq.e) 
  
# write.csv(res, file = "output/sim2_cjc.csv")


  
# END
