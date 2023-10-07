# Oct 4, 2023
# Started by D Loughnan
#Aim of this model --- to check whether model with lambda out performs model without
# This is the lambda model--- lambda on the slope and intercept---one cue (x)
# Testing different values lambda (lambda slope = lambda intercept): 0, 0.2,0.8

if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("Lizzie", getwd())>0)) {
    setwd("~/Documents/git/teaching/stan/pmm/sandbox")
} else{
    setwd("/home/deirdre/pmm") # for midge
}

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)

rm(list=ls())
# Set number of cores available
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

### Simulate data with and without lambda
nspecies = 5
nind = 5

ndat <- nspecies*nind

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
nruns <- 2
nparamL <- 8 # number of parameters including lp__
ntotL <- nparamL+2*nspecies

nparamN <- 6 # no lambda parameters
ntotN <- nparamN+2*nspecies
lambda = c(0, 0.2 ,0.8)
#lambda = c(0,0,0.8,0.8)
nlambda <- length(lambda)

sumDatL <- data.frame(matrix(NA, ntotL*nruns*nlambda, 1))
names(sumDatL) <- c("rep")
sumDatL$rep <- rep(1:ntotL, nruns*nlambda )
sumDatL$runs <- rep(1:nruns, each = ntotL)
sumDatL$lambda <- rep(lambda, each = ntotL*nruns)

sumDatN <- data.frame(matrix(NA, ntotN*nruns*nlambda, 1))
names(sumDatN) <- c("rep")
sumDatN$rep <- rep(1:ntotN, nruns*nlambda )
sumDatN$runs <- rep(1:nruns, each = ntotN)
sumDatN$lambda <- rep(lambda, each = ntotN*nruns)

sumDatLambOut <- vector()

sumDatNoOut <- vector()

for (l in 1:length(lambda)){
  
  for (j in 1:nruns){
param <- list(a_z = 4, # root value intercept
              lam_interceptsa = lambda[l], # lambda intercept
              sigma_interceptsa = 0.2, # rate of evolution intercept
              b_zf = 0.6, # root value trait1 slope
              lam_interceptsbf = lambda[l], # lambda trait1
              sigma_interceptsbf = 0.1, # rate of evolution trait1
              sigma_y = 0.01 # overall sigma
              )
 param

# Generate data for zero lambda data set
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsbf"]])         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)

dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind),
                       intercept = rep(intercepts[i], nind),
                       x1 = rnorm(n = nind, mean = 10, sd = 3),
                       trait1 = rep(slopes_bf[i], nind))
    dfhere <- rbind(dfhere, temp)
}
dfhere$mu <- dfhere$intercept + dfhere$x1 * dfhere$trait1
dfhere$y <- rnorm(n = nrow(dfhere), mean = dfhere$mu, sd = param[["sigma_y"]])

write.csv(dfhere, paste("output/simData", lambda[l], j, ".csv", sep =""))

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
    a_z.temp <- rnorm(n = nspecies, mean = param[["a_z"]], sd = 1)
    b_z.temp <- rnorm(n = nspecies, mean = param[["b_zf"]], sd = 1)
    return(append(list(a = a_z.temp, b_force = b_z.temp),
                  param))
}

mdl.data <- list(y = dfhere$y,
                 N = nrow(dfhere),
                 n_sp = nspecies,
                 sp =dfhere$sp,
                 x1=dfhere$x1,
                 Vphy=vcv(spetree, corr = TRUE))
                 
testLamb <- stan("stan/phyloMdlLambdaIntSlope.stan", 
                 #  control = list(max_treedepth =15),
                 data = mdl.data,
                 #init = simu_inits,
                 iter = 4000,
                 warmup = 3000,
                 chains = 4)

save(testLamb, file = paste("output/testLambdaEffect_", lambda[l], j,".Rda", sep =""))

#sumerL <- summary(testLamb)$summary

sumerLdf <- data.frame(summary(testLamb)$summary)
sumDatLambOut <- rbind(sumDatLambOut, sumerLdf)



testNoLamb <- stan("stan/phyloMdlNoLambda.stan",
                 #  control = list(max_treedepth =15),
                 data = mdl.data,
                 #init = simu_inits,
                 iter = 4000,
                 warmup = 3000,
                 chains = 4)
save(testNoLamb, file = paste("output/testNoLambdaEffect_", lambda[l], j,".Rda", sep =""))
sumerN <- summary(testNoLamb)$summary

sumerNdf <- data.frame(summary(testNoLamb)$summary)
sumDatNoOut <- rbind(sumDatNoOut, sumerNdf)
#muTraitSp <- sumer[grep("a\\[", rownames(sumer))]
# pdf("LambdaSpComp0.2.pdf", height = 5, width = 5)
# plot(muTraitSp ~ intercepts, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# dev.off()


  }
}
write.csv(sumDatLambOut, paste("output/mdlOutLambdaReppedTemp.csv", sep = ""))
write.csv(sumDatNoOut, paste("output/mdlOutNoLambdaReppedTemp.csv", sep = ""))

sumDatL <- cbind(sumDatL, sumDatLambOut)
write.csv(sumDat, paste("output/mdlOutLambdaRepped.csv", sep = ""))

sumDatN <- cbind(sumDatN, sumDatNoOut)
write.csv(sumDatN, paste("output/mdlOutNoLambdaRepped.csv", sep = ""))
# sumer <- data.frame(summary(test_new)$summary[c("a_z","lam_interceptsa","sigma_interceptsa", "b_z","lam_interceptsb" ,"sigma_interceptsb","sigma_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
# sumer <- data.frame(sumer)
# sumer$param <- do.call(rbind.data.frame, t(param))
# colnames(sumer)[colnames(sumer) == "c.4..0..0.2..0.6..0..0.1..0.01."] <- "testValue"

##########################################################################################
# Compare the above to a model without any lambda
# intercepts <- rnorm(nspecies, param[["a_z"]], param[["sigma_interceptsa"]] )
# # Generate bf slope
# slopes_bf <- rnorm(nspecies, param[["b_zf"]],  param[["sigma_interceptsbf"]])
# 
# dfNoLam <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
# for (i in 1:nspecies){
#     temp <- data.frame(sp = rep(i, nind),
#                        intercept = rep(intercepts[i], nind),
#                        x1 = rnorm(n = nind, mean = 10, sd = 3),
#                        trait1 = rep(slopes_bf[i], nind))
#     dfNoLam <- rbind(dfNoLam, temp)
# }
# dfNoLam$mu <- dfNoLam$intercept + dfNoLam$x1 * dfNoLam$trait1
# dfNoLam$y <- rnorm(n = nrow(dfNoLam), mean = dfNoLam$mu, sd = param[["sigma_y"]])
# 
# 
# 
# mdl.data <- list(yobs = dfNoLam$y,
#                  N = nrow(dfNoLam),
#                  Nspp = nspecies,
#                  species =dfNoLam$sp,
#                  x1= dfNoLam$x1
#                  )
# 
# test_new <- stan("sandbox/stan/phyloMdlNoLambda.stan", 
#                  #  control = list(max_treedepth =15),
#                  data = mdl.data,
#                  #init = simu_inits,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4)
# 
# sumer <- summary(test_new)$summary
# save(test_new, file = paste("output/testLambdaEffect_noLambda.Rda", sep =""))
# 
# 
# muTraitSp <- sumer[grep("a_sp\\[", rownames(sumer))]
# pdf("noLambdaSpComp.pdf", height = 5, width = 5)
# plot(muTraitSp ~intercepts, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# dev.off()
# 
# sumerdf <- data.frame(summary(test_new)$summary)
# write.csv(sumerdf, paste("output/mdlOutNoLambda.csv"))
# 
# nolam <- read.csv("output/mdlOutNoLambda.csv")
# names(nolam) <- c("parameters","mean","se_mean", "sd","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")
# 
# paramNames <- c("mu_a", "mu_b","sigma_a","sigma_b","sigma_y")
# values <- c("mean","2.5%","25%","50%", "75%","97.5%")
# 
# nolamSum <- nolam[nolam$parameters %in% paramNames, c("parameters","mean","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")]
# nolamSum$testV <- c(param[["sigma_y"]],param[["a_z"]],param[["sigma_interceptsa"]],param[["b_zf"]],param[["sigma_interceptsbf"]])
# nolamSum
# 
# # sumer <- data.frame(summary(test_new)$summary[c("a_z","lam_interceptsa","sigma_interceptsa", "b_z","lam_interceptsb" ,"sigma_interceptsb","sigma_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
# # sumer <- data.frame(sumer)
# # sumer$param <- do.call(rbind.data.frame, t(param))
# # colnames(sumer)[colnames(sumer) == "c.4..0..0.2..0.6..0..0.1..0.01."] <- "testValue"
# # write.table(sumer, file = paste("output/testLambdaEffect_", param$lam_interceptsa, ".csv"))
# # 
# # save(test_new, file = paste("output/testLambdaEffect_", param$lam_interceptsa, ".Rda"))
# 
# #######################################################################################
# 
# # parameters
# param <- list(a_z = 4, # root value intercept
#               lam_interceptsa = 0, # lambda intercept
#               sigma_interceptsa = 0.2, # rate of evolution intercept
#               b_zf = 0.6, # root value trait1 slope
#               lam_interceptsbf = 0, # lambda trait1
#               sigma_interceptsbf = 0.1, # rate of evolution trait1
#               sigma_y = 0.01 # overall sigma
# )
# 
# 
# 
# lam0 <- read.csv("output/mdlOutLambda0.csv")
# #head(lam0)
# names(lam0) <- c("parameters","mean","se_mean", "sd","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")
# 
# paramNames <- c("a_z","lam_interceptsa","sigma_interceptsa", "b_z","lam_interceptsb" ,"sigma_interceptsb","sigma_y")
# values <- c("mean","2.5%","25%","50%", "75%","97.5%")
# 
# lam0Sum <- lam0[lam0$parameters %in% paramNames, c("parameters","mean","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")]
# lam0Sum$testV <- c(param[["sigma_y"]],param[["lam_interceptsa"]],param[["sigma_interceptsa"]],param[["lam_interceptsbf"]],param[["sigma_interceptsbf"]],param[["b_zf"]],param[["a_z"]])
# lam0Sum
# ### Lambda = 0.2
# 
# param <- list(a_z = 4, # root value intercept
#               lam_interceptsa = 0.2, # lambda intercept
#               sigma_interceptsa = 0.2, # rate of evolution intercept
#               b_zf = 0.6, # root value trait1 slope
#               lam_interceptsbf = 0.2, # lambda trait1
#               sigma_interceptsbf = 0.1, # rate of evolution trait1
#               sigma_y = 0.01 # overall sigma
# )
# 
# lam0.2 <- read.csv("output/mdlOutLambda0.2.csv")
# names(lam0.2) <- c("parameters","mean","se_mean", "sd","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")
# 
# #head(lam0.2)
# lam0.2Sum <- lam0.2[lam0.2$parameters %in% paramNames, c("parameters","mean","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")]
# lam0.2Sum$testV <- c(param[["sigma_y"]],param[["lam_interceptsa"]],param[["sigma_interceptsa"]],param[["lam_interceptsbf"]],param[["sigma_interceptsbf"]],param[["b_zf"]],param[["a_z"]])
# lam0.2Sum
# ### Lambda = 0.8
# 
# param <- list(a_z = 4, # root value intercept
#               lam_interceptsa = 0.8, # lambda intercept
#               sigma_interceptsa = 0.2, # rate of evolution intercept
#               b_zf = 0.6, # root value trait1 slope
#               lam_interceptsbf = 0.8, # lambda trait1
#               sigma_interceptsbf = 0.1, # rate of evolution trait1
#               sigma_y = 0.01 # overall sigma
# )
# 
# lam0.8 <- read.csv("output/mdlOutLambda0.8.csv")
# names(lam0.8) <- c("parameters","mean","se_mean", "sd","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")
# 
# lam0.8Sum <- lam0.8[lam0.8$parameters %in% paramNames, c("parameters","mean","X2.5.","X25.","X50.", "X75.","X97.5.","n_eff","Rhat")]
# lam0.8Sum$testV <- c(param[["sigma_y"]],param[["lam_interceptsa"]],param[["sigma_interceptsa"]],param[["lam_interceptsbf"]],param[["sigma_interceptsbf"]],param[["b_zf"]],param[["a_z"]])
# lam0.8Sum





      
