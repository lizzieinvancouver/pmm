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

nspecies = 40
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# now set up the trait
m <- 0.6
lam <- 0.8
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

if(FALSE){ # alternative test data (NOT DONE)
Vphy <- vcv(spetree)
diag(Vphy) <- 0

slopez <- MVN(0, Vphy) # START here! (and be sure to review all of the below as I just pasted it in for now)
# rep_vector(0,n_sp), diag_matrix(rep_vector(null_interceptsb, n_sp)) + lam_interceptsb*Vphy

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

}

# The below two lines make a big improvement (compared to Vphy=vcv(spetree, corr=TRUE))
# but ... they seem to just divide the lambda that used to be almost always 1 to alsmost always 0.5 now
# always 0.5 now as lam.int and null.int are always the same value,
# but when Vphy=vcv(spetree, corr=TRUE), then lam.int >> null.int and estimate approaches 1
Vphy=vcv(spetree)
diag(Vphy) <- 0

testme <- stan("stan/ubermini.stan",
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=Vphy), # Vphy=vcv(spetree, corr=TRUE)),
                iter=2000, chains=4, seed=123456)
summary(testme)$summary

# From Will's code
# ...the estimated phylogenetic signal (Pagel's Lambda) of the
# estimated environmental response is ~.7 (i.e., there's signal)
# (see Hadfield & Nakagawa (2010) J. Evol. Biol. 23(3): 494-508
# (bottom of pg 497, which is also just the classic from Housworth etc.) 
# for this method of calculating it)
lam.int <- mean(extract(testme)[["lam_interceptsb"]])
null.int <- mean(extract(testme)[["null_interceptsb"]])
lam.int / (null.int + lam.int) 

sumer <- summary(testme)$summary
# sumer[grep("b_force", rownames(sumer)),]

# Compare true slopes to estimated slopes
sumer[grep("b_force", rownames(sumer)), "mean"]
slopez

plot(slopez, sumer[grep("b_force", rownames(sumer)), "mean"])
