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

nspecies = 90
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait
m <- 0.6
lam <- 0.7
sigy <- 0.01
sig2 <- 0.1

scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sig2)
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
    
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigy)

testme <- stan("stan/ubermini_2.stan", # Note: changed to a new model
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=vcv(spetree)),# Note: dropped the corr=TRUE (is that correct? Lizzie says: Doesn't seem to matter, sadly)
                iter=2000, chains=2, seed=123456)
summary(testme)$summary
summary(testme)$summary[c("lam_interceptsb","sigma_interceptsb","b_z"),"50%"]
fitContinuous(spetree, slopez, model="lambda")

# Compare true slopes to estimated slopes
sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)), "mean"]
# slopez

plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])
