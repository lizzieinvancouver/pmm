## Started 3 Sep 2020 ##
## Generate test data for my phylogeny-on-the-slope models ##

# y ~ b_force[phylo]*x1 + error

# Setting working directory
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")

## load packages
require(phytools)
library(MASS)
require(rstan)

nspecies = 10
nind = 3

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE, scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait and the x (co-factor) basics
m <- 0.6
sigy <- 0.01
xmean <- 10
xsd <- 3

# Create the slopes without any fancy phy packages
Vphy <- vcv(spetree, corr=TRUE)
null_interceptsb <- 0.5
lam_interceptsb <- 0.01

slopez <- mvrnorm(n=1, mu=rep(0, nspecies), Sigma=(diag(rep(null_interceptsb, nspecies)) + lam_interceptsb*Vphy))

# Now, we have the slopes, so build the rest of the data
dfhere <- data.frame(x=numeric(), y=numeric())

for (i in 1:length(slopez)){
    slopehere <- slopez[i]
    xhere <- rnorm(nind, xmean, xsd) # these are experiments, no phylo structure in x 
    yhere <- xhere*slopehere
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}

dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
    
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigy)

# Run the model
testme <- stan("stan/ubermini.stan",
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=vcv(spetree, corr=TRUE)),
                iter=5000, chains=4, seed=123456)
summary(testme)$summary

sumer <- summary(testme)$summary

# Compare true slopes to estimated slopes
plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])

# Compare matrix parameters
lam.int <- mean(extract(testme)[["lam_interceptsb"]])
null.int <- mean(extract(testme)[["null_interceptsb"]])
lam.int
null.int
# true values were ....
null_interceptsb
lam_interceptsb

