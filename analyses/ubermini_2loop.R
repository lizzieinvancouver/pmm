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

options(mc.cores = parallel::detectCores())

dfout <- data.frame(m=numeric(), lam=numeric(), phylosig=numeric(),
    bz50=numeric(), bz2.5=numeric(), bz97.5=numeric(),
    lam_interceptsb50=numeric(), lam_interceptsb2.5=numeric(), lam_interceptsb97.5=numeric(),
    slopelm=numeric(), Rhat=numeric())

# Set up one species and tree for all sims ... 
nspecies = 90
nind = 10

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

msims <- c(0.1, 0.5, 1, 5)
lamsims <- c(0.1, 0.5, 0.9)

for(i in c(msims)){
    for (j in c(lamsims)){

# now set up the trait
m <- i
lam <- j
sigy <- 0.01
sig2 <- 0.1

scaledtree <- rescale(spetree, model="lambda", lam)
slopez <- fastBM(scaledtree, a=m, mu=0, sig2=sig2)
phylosig(x=slopez, tree=spetree, method="lambda")

# for testing ...
nulltree <- rescale(spetree, model="lambda", 0)

dfhere <- data.frame(x=numeric(), y=numeric())

for (k in 1:length(slopez)){
    slopehere <- slopez[k]
    xhere <- rnorm(nind, 10, 3) # these are experiments, no phylo structure in x 
    yhere <- xhere*slopehere
    dfadd <- data.frame(x=xhere, y=yhere)
    dfhere <- rbind(dfhere, dfadd)
}

dfhere$sp <- rep(spetree$tip.label, each=nind)
dfhere$spnum <- as.numeric(gsub('s', '', dfhere$sp))
    
dfhere$yerr <- dfhere$y + rnorm(nrow(dfhere), 0, sigy)

testme <- stan("stan/ubermini_2_biggerpriors.stan", # Note: changed to a new model
                data=list(N=nrow(dfhere), n_sp=nspecies, sp=dfhere$spnum,
                x=dfhere$x, y=dfhere$yerr,
                Vphy=vcv(spetree)), # Note: dropped the corr=TRUE
                iter=4000, chains=4)


# Compare true slopes to estimated slopes
sumer <- summary(testme)$summary
        
slopemodel <- summary(lm(slopez ~ sumer[grep("b_force", rownames(sumer)), "mean"]))

dfadd <- data.frame(m=m, lam=lam, phylosig=phylosig(x=slopez, tree=spetree, method="lambda")[[1]],
    bz50=summary(testme)$summary[c("b_z"),"50%"],
    bz2.5=summary(testme)$summary[c("b_z"),"2.5%"],
    bz97.5=summary(testme)$summary[c("b_z"),"97.5%"],
    lam_interceptsb50=summary(testme)$summary[c("lam_interceptsb"),"50%"],
    lam_interceptsb2.5=summary(testme)$summary[c("lam_interceptsb"),"97.5%"],
    lam_interceptsb97.5=summary(testme)$summary[c("lam_interceptsb"),"2.5%"],
    slopelm=coef(slopemodel)[[2]],
    Rhat=mean(summary(testme)$summary[c("lam_interceptsb", "sigma_interceptsb", "b_z", "sigma_y"),"Rhat"]))

dfout <- rbind(dfout, dfadd)

    }
}

write.csv(dfout, "output/dfout_nocorr.csv")

# df <- read.csv("~/Documents/git/teaching/stan/pmm/analyses/output/dfout_nocorr.csv")
# ggplot(df, aes(phylosig, lam_interceptsb50, group=as.factor(m), color=as.factor(m)))+geom_point()
# ggplot(df, aes(m, bz50, group=as.factor(m), color=as.factor(lam)))+geom_point()
