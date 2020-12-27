## Started 28 November 2020 ##
## By Lizzie ## 
## Taken mainly from Phylo_ospree_reanalyses.R ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/teaching/stan/pmm/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

# Get the data
# See Phylo_ospree_reanalyses.R for what I wrote out  (ospree repo)

library("ape")
library("rstan")

phylo <- read.tree("input/ospreephylo.tre")
bb.stan <- read.csv("input/ospreebb.csv")

# Step 1: Get spps and VCVPHY in same order
# bb.stan$spps[phylo$tip.label]
phylo$tip.label
d <- bb.stan[match(phylo$tip.label, bb.stan$spps),]

phymatch <- data.frame(tip=phylo$tip.label, sppnum=c(1:length(phylo$tip.label)))
d <- merge(bb.stan, phymatch, by.x="spps", by.y="tip")
d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)

if(FALSE){
# the below runs ...
testme <- stan("stan/ubermini_2_biggerpriors.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x=d$force.z, # chill=d$chill.z, photo=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=3000, chains=4)

# ran this with one chain, 4 div trains and some structure to sigma_interceptsb
testme2 <- stan("stan/lessmini_oneslopeintercept.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x=d$force.z, # chill=d$chill.z, photo=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=3000, chains=1)


# ran this with one chain, 39 div trains and some structure to sigma_interceptsbf and lam_interceptsbf
testme3 <- stan("stan/lessmini_twoslopeintercept.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x1=d$force.z, x2=d$chill.z, # photo=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=3000, chains=1)
}

# 312 divergent transitions, some structure in sigma_interceptsbf and sigma_interceptsbp
testme4 <- stan("stan/lessmini_threeslopeintercept.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x1=d$force.z, x2=d$chill.z, x3=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=3000, chains=3)


if(FALSE){
## ...
testme <- stan("stan/nointer_2levelphyall_2.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       force=d$force.z, chill=d$chill.z, photo=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=3000, chains=1)
}

# goo <- readRDS("output/testme4.Rds")
# on server: 3500/5000: 338 DivTrans 5500/7000: 204 DivTrans; 5500/7000 and 0.9: 134
# testme4 (new) 8000/9000 -- 227; testme 5: 9000/10000 and 0.9 -- 137; testme6 9000/10000 and 0.95 -- 197; testme7: 9000/10000 and 0.99 -- 42!
# next on server (took at least 2 days I think!): testme4 8000/9000 and 0.99 -- 32 DivTrans; testme5  11000/12000 and 0.99 24; testme6 8000/9000 and 0.999 37;  testme 7 11000/12000 and 0.999 14 DivTrans. 
goo <- readRDS("output/testme4.Rds")

