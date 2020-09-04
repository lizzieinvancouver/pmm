## Started 3 Sep 2020 ##
## Generate test data for my phylogeny-on-the-slope models ##

# y ~ b_force[phylo]*x1 + error

# Setting working directory
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")

## load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

nspecies = 10
nind = 3

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait and the x (co-factor) basics
m <- 0.6
sigy <- 0.01
xmean <- 10
xsd <- 3

# set up adjusting the trait by lambda in geiger
lam <- 0.7
sig2 <- 0.1
    
scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sig2)
phylosig(x=slopez, tree=spetree, method="lambda")

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

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
phylosig(x=slopez, tree=spetree, method="lambda")

null.int # in this version of the code, this parameter is not clearly set (because I set up lambda in geiger)
