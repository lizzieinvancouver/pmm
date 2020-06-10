## Started 9 June 2020 ##
## Generate test data for my phylogeny Stan models ##

## Someday we want:
# y ~ a[sp] + b_force[phylo]*x1 + b_chill[phylo]*x2 + b_photo[phylo]*x3 + error

## Almost simple to start ...
## skip an intercept and just estimate
# y ~ a[sp] + b_force*x1 + error

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

nspecies = 40
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the intercepts
a <- 10
m <- 0.6
lam <- 0.5
sigy <- 0.01
sig2 <- 0.1

scaledtree <- rescale(spetree, model="lambda", lam)
interceptz <- fastBM(scaledtree, a=a, mu=0, sig2=sig2)
phylosig(x=interceptz, tree=spetree, method="lambda")


# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

dfhere <- data.frame(x=numeric(), y=numeric())

for (i in 1:length(interceptz)){
    intercepthere <- interceptz[i]
    xhere <- rnorm(nind, 10, 3) # these are experiments, no phylo structure in x 
    yhere <- intercepthere + xhere*m
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}

dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
    
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigy)

testme <- stan("stan/uberlessmini.stan",
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=vcv(spetree, corr=TRUE)),
                iter=2000, chains=4, seed=123456)
summary(testme)$summary

# From Will's code
# ...the estimated phylogenetic signal (Pagel's Lambda) of the
# estimated environmental response is ~.7 (i.e., there's signal)
# (see Hadfield & Nakagawa (2010) J. Evol. Biol. 23(3): 494-508 for
# this method of calculating it)
lam.int <- mean(extract(testme)[["lam_interceptsb"]])
null.int <- mean(extract(testme)[["null_interceptsb"]])
lam.int / (null.int + lam.int) 

sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)),]

# Compare true slopes to estimated slopes
sumer[grep("b_force", rownames(sumer)), "mean"]
interceptz

plot(interceptz, sumer[grep("a\\[", rownames(sumer)), "mean"])
