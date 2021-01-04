## Started 28 November 2020 ##
## By Lizzie ## 
## Taken mainly from Phylo_ospree_reanalyses.R ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

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

# Step 2: Run the model
testme6 <- stan("stan/lessmini_threeslopeintercept.stan",
    data=list(N=nrow(d), n_sp=nspecies, sp=d$sppnum,
       x1=d$force.z, x2=d$chill.z, x3=d$photo.z,
       y=d$resp, Vphy=vcv(phylo)),
    iter=9000, warmup=8000, chains=4, cores=4, control=list(adapt_delta=0.999))

saveRDS(testme6, "output/testme6.rds")
