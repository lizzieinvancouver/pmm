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
require(phyr)

nspecies = 90
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
slopez <- fastBM(scaledtree, a=0, mu=0, sig2=sig2)
phylosig(x=slopez, tree=spetree, method="lambda")


## quick check for how BM traits are simulated
sigs<-list()
for(i in 1:100){
  print(i)
  slopez <- fastBM(scaledtree, a=0, mu=0, sig2=sig2)
  sigs[[i]]<-phylosig(x=slopez, tree=spetree, method="lambda")$lambda
}
hist(unlist(sigs))
mean(unlist(sigs))
  
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


## add species names and errors to simulated data
dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigy)


# check whether simulated x (should not have phylo signal) actually does not have it
checkx <- aggregate(dfhere, by=list(sp=dfhere$sp),mean)
phylosig(x=checkx$x, tree=spetree, method="lambda")$lambda



# Run the model
# see https://cran.r-project.org/web/packages/phyr/phyr.pdf pg 23-24
Vphy=vcv(spetree)


## first model specification is the simplest possible 
## partial pooling on slopes for the phylogeny with no fixed effects
testme <- pglmm(y ~ x|sp, data = dfhere,
                family = "gaussian", cov_ranef = list(sp = Vphy))
testme

## The above specification does the partial pooling by species, disregarding
## the phylogeny. This can be checked by inspecting the VCV matrix used for
## the random effects (which is an identity matrix, 0s and 1s in the diags).  
image(testme$random.effects$`x|sp`$covar)


## to specify the above model to actually account for the phylogeny, 
## underscores must be used
testme.phy <- pglmm(y ~ x|sp__, data = dfhere,
                    family = "gaussian", 
                    cov_ranef = list(sp = Vphy))

## now there are two random effects |sp (species, non phylogenetic) and
## |sp__ (phylogenetic) 
image(testme.phy$random.effects$`x|sp`$covar) # same as before
image(testme.phy$random.effects$`x|sp__`$covar) # phylogenetic structure


## now we would ideally test lambda or some sort of phylogenetic signal (a cool
## thing from lamda is that we can compare directly against other models)
## following the logic of the BRMS vignette (by Paul Burkner) an analogous 
## metric of phylogenetic signal would be calcultated as the ratio between 
## variance attributable to the phylogeny alone and total variance: 

phylo.var <- testme.phy$s2r[2]/sum(testme.phy$s2r,testme.phy$s2resid)


## would it be possible to obtain a 'lambda' metric (multiplier of VCV off-diagonal
## elements) that is comparable to that in other models (pgls, stan gaussian, etc.)


