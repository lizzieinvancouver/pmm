## Started 29 May 2020 ##
## Generate test data for my phylogeny Stan models ##

## Someday we want:
# y ~ a[sp] + b_force[phylo]*x1 + b_chill[phylo]*x2 + b_photo[phylo]*x3 + error

## Uber simple to start ...
## skip an intercept and just estimate
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

# Set one of the below to TRUE
geigerversion=FALSE # creates test data using Brownian motion model the phylo-folks like
mvnversion=TRUE # creates test data to match Stan code (I think)

nspecies = 40
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait
m <- 0.6
sigy <- 0.01
xmean <- 10
xsd <- 3

if(geigerversion){
lam <- 1
sig2 <- 0.1
    
scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sig2)
phylosig(x=slopez, tree=spetree, method="lambda")

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

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
}

if(mvnversion){ 
library(MASS)

Vphy <- vcv(spetree)
diag(Vphy) <- 0
null_interceptsb <- 0.5
lam_interceptsb <- 0.01

slopez <- mvrnorm(n=1, mu=rep(0, nspecies), Sigma=(diag(rep(null_interceptsb, nspecies)) + lam_interceptsb*Vphy))

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

}

# Non-identifiable? Using geigerversion and the below two lines compared to Vphy=vcv(spetree, corr=TRUE) ...
# seems to just divide the lambda that used to be almost always 1 to alsmost always 0.5 now ...
# as lam.int and null.int are always the same value,
# but when Vphy=vcv(spetree, corr=TRUE), then lam.int >> null.int and estimate approaches 1
Vphy=vcv(spetree)
diag(Vphy) <- 0

testme <- stan("stan/ubermini.stan",
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=Vphy), # Vphy=vcv(spetree, corr=TRUE)),
                iter=5000, chains=4, seed=123456)
summary(testme)$summary

sumer <- summary(testme)$summary

# Compare true slopes to estimated slopes
sumer[grep("b_force", rownames(sumer)), "mean"]
slopez

plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])

lam.int <- mean(extract(testme)[["lam_interceptsb"]])
null.int <- mean(extract(testme)[["null_interceptsb"]])

lam.int
null.int

# Will suggests ...
lam.int / (null.int + lam.int) 


# From Will's code
# ...the estimated phylogenetic signal (Pagel's Lambda) of the
# estimated environmental response is ~.7 (i.e., there's signal)
# (see Hadfield & Nakagawa (2010) J. Evol. Biol. 23(3): 494-508
# (bottom of pg 497, which is also just the classic from Housworth etc.) 
# for this method of calculating it)
lam.int / (null.int + lam.int) 


if(FALSE){
# From Geoff on 31 May 2020:
# Here's one set of calls that might clarify what's happening:

# Regular diagonal
diag(rep(1, 10))

# Add a constant
diag(rep(1, 10)) + 6

# Add a vector
diag(rep(1, 10)) + c(1:10)
}
