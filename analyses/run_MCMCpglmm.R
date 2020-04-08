## Code from Joly et al. 2019 ##
## He sent this in an email on 10 March 2020 ##

## Lizzie asked: Do you have any simulated data + model code from the MEE paper that using just an interspecific phylogeny? Then we could do some work to better understand what the within-species centering is doing.

## He wrote: I don't have code for exactly what you want, but you can probably adapt the one I have relatively easily. In the R code attached, the function "one.sim.pmm()" simulates data and analyze it with PMM and PGLS. I think I used the same notation as in the paper. You would just have to drop the term for the intraspecific variance and simulate more values for the residual error (or add a similar term for the measurement error in addition to the residual error). Let me know if you have any other questions.

## Modified by Nacho on 8th april - to source functions to simulate data and run MCMCglmm and PGLS models
## 

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} 

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



## Parameters (set your own set of perameters according to PMM_models_comparison.csv)
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

simdata<-one.sim.pmm(nspecies = nspecies, nindividuals = nindividuals,B = B,
                     sigma.sq.x = sigma.sq.x,sigma.sq.p = sigma.sq.p,
                     sigma.sq.c = sigma.sq.c,sigma.sq.e = sigma.sq.e)

  
#plot(simdata$data$x,simdata$data$y)  
 

## run MCMCglmm models 
# Prepare infraspecific correlation matrix for MCMCglmm
#
# Single value decomposition of the intraspecific structure
#   -> follow code of Stone et al. 2012
    intra.svd <- svd(simdata$intramat)
    intra.svd <- intra.svd$v %*% (t(intra.svd$u) * sqrt(intra.svd$d))
    rownames(intra.svd) <- colnames(intra.svd) <- rownames(simdata$intramat)

# The following is important and MCMCglmm searches the variables in the 
# global environment
    intra.svd <<- intra.svd
    
###
# Run the MCMCpglmm analyses with different correlation structures.
# For each, we will use a very diffuse prior.
    
# Model M.0 -> no random effects
    priorpr.m0 <- list(R = list(V = 1, nu = 0.002))
    M.0 <- MCMCglmm(y ~ x, 
                    data=simdata$data, scale=TRUE,prior=priorpr.m0,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
    
# Model M.1 -> only phylogenetic random effects
    priorpr.m1 <- list(R = list(V = 1, nu = 0.002), 
                       G = list(G1 = list(V = 1, nu = 0.002)))
    M.1 <- MCMCglmm(y ~ x, random = ~ animal, pedigree = simdata$phylo,
                    data=simdata$data, scale=TRUE,prior=priorpr.m1,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
  
    
## run PGLS models 
#

inter.mat.inflated <- inflate.mat(simdata$intermat,nspecies,nindividuals)

# Generate phylogenetic correlation matrices
# The intra correlation structure takes the inter-specific and the intra-
# specific correlation matrices as input.
pgls.intra.cor <- corIntra(0, inter=inter.mat.inflated, intra=simdata$intramat, fixed=TRUE)
pgls.inter.cor <- corIntra(1, inter=inter.mat.inflated, intra=simdata$intramat, fixed=TRUE)
pgls.intrainter.cor <- corIntra(0.5, inter=inter.mat.inflated, intra=simdata$intramat, fixed=FALSE)
    
# Fit PGLS
gls.mod0 <- gls.mod1 <- gls.mod2 <- gls.mod3 <-NULL
try(gls.mod0 <- gls(y~x,data=simdata$data)) # no phylogenetic correlation
try(gls.mod1 <- gls(y~x,correlation = pgls.inter.cor, data=simdata$data),silent=TRUE) # inter structure 
try(gls.mod2 <- gls(y~x,correlation = pgls.intra.cor, data=simdata$data),silent=TRUE) # intra structure
try(gls.mod3 <- gls(y~x,correlation = pgls.intrainter.cor, data=simdata$data),silent=TRUE) # Estimating delta 
    

###
# Compile the results
  
# Model summaries
  M.0.sum <- summary(M.0)
  M.1.sum <- summary(M.1)
#  M.2.sum <- summary(M.2)
#  M.3.sum <- summary(M.3)
  gls.mod0.sum <- summary(gls.mod0)
  gls.mod1.sum <- summary(gls.mod1)
  gls.mod2.sum <- summary(gls.mod2)
  gls.mod3.sum <- summary(gls.mod3)
  
# Data frame
# 
  
  res <- data.frame(models=c("M.0","M.1","gls0","gls1","gls2","gls3"),
                    #random.effects=c("NA","inter","infra","inter+infra","ols","glsInter","glsIntra","glsInterIntra"),
                    slope=c(M.0.sum$solutions[2,1],M.1.sum$solutions[2,1],
                            gls.mod0$coefficients[2],gls.mod1$coefficients[2],
                            gls.mod2$coefficients[2],gls.mod3$coefficients[2]),
                    slope.sd=c(sd(M.0$Sol[,'x']),sd(M.1$Sol[,'x']),
                               gls.mod0.sum$tTable[2,2],
                               gls.mod1.sum$tTable[2,2],
                               gls.mod2.sum$tTable[2,2],
                               gls.mod3.sum$tTable[2,2]),
                    slope.signif=c(ifelse(M.0.sum$solution[2,5]<0.05,1,0),
                                   ifelse(M.1.sum$solution[2,5]<0.05,1,0),
                                   ifelse(gls.mod0.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod1.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod2.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod3.sum$tTable[2,4]<0.05,1,0)),
                    DIC=c(M.0$DIC,M.1$DIC,-999,-999,-999,-999))
  
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
  
  
  # Heritability (h^2)
  res$h2 <- c(NA,M.1.var["animal"]/sum(M.1.var),
              NA,NA,NA,NA)
  
#write.csv(res,file = "output/sim1_IMC.csv")

  
# END
