## Started 5 July 2020 ##
## Generate test pgls ##


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
nind = 1

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# attempt to create the model with phylogenetic structure in the error
a <- 10
m <- 0.6
lam <- 0.5
sigy <- 0.1
sig2 <- 0.1

# we want the structure on the error, hmmm....
scaledtree <- rescale(spetree, model="lambda", lam)
errorz <- fastBM(scaledtree, a=sigy, mu=0, sig2=sig2)
phylosig(x=errorz, tree=spetree, method="lambda")

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

dfhere <- data.frame(x=numeric(), y=numeric())

for (i in 1:nspecies*nind){
    xhere <- rnorm(nind, 10, 3) # these are experiments, no phylo structure in x 
    yhere <- a + xhere*m
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}

dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))

# Not at all sure this is correct ...
dfhere$yerr <- dfhere$y + errorz

# Construct Lmat
Lmat <- matrix(rep(1), nrow = nspecies, ncol = nspecies)
diag(Lmat) <- 0


# Sadly, not running because it's "Rejecting initial value" ... no time to fix just now
testme <- stan("stan/uberpgls.stan",
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=vcv(spetree),
                Lmat=Lmat),
                iter=2000, chains=4)
summary(testme)$summary


sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)),]
