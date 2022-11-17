## Started 9 September 2021 ##
## By Lizzie ##

## Comparing my model (Jan 2021) using OSPREE data to Geoff's version (July 2021) 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan)

setwd("~/Documents/git/teaching/stan/pmm/analyses/")

oldmod <- readRDS("output/testme8.rds")
newmod <- readRDS("output/testme.rds")

source("source/stan_utility.R")
check_all_diagnostics(oldmod)
check_all_diagnostics(newmod)

old <- extract(oldmod)
new <- extract(newmod)

oldsum <- summary(oldmod)$summary
newsum <- summary(newmod)$summary

aold <- oldsum[grep("a\\[", rownames(oldsum)),]
anew <- newsum[grep("a\\[", rownames(newsum)),]

par(mfrow=c(2,2))
plot(aold[,1]~anew[,1], ylab=" intercepts old model", xlab="intercepts new model")
abline(0, 1, col="dodgerblue")
plot(oldsum[grep("b\\_photo\\[", rownames(oldsum)),][,1]~newsum[grep("b\\_photo\\[", rownames(newsum)),][,1],
     ylab="photo slope old model", xlab="photo slope new model")
abline(0, 1, col="dodgerblue")
plot(oldsum[grep("b\\_force\\[", rownames(oldsum)),][,1]~newsum[grep("b\\_force\\[", rownames(newsum)),][,1],
     ylab="force slope old model", xlab="force slope new model")
abline(0, 1, col="dodgerblue")
plot(oldsum[grep("b\\_chill\\[", rownames(oldsum)),][,1]~newsum[grep("b\\_chill\\[", rownames(newsum)),][,1],
     ylab="chill slope old model", xlab="chill slope new model")
abline(0, 1, col="dodgerblue")

oldsum[grep("lam", rownames(oldsum)),]
newsum[grep("lam", rownames(newsum)),]

oldsum[grep("sigma", rownames(oldsum)),]
newsum[grep("sigma", rownames(newsum)),]
