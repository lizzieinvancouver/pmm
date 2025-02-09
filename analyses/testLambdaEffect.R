# Oct 4, 2023
# Started by D Loughnan
#Aim of this model --- to check whether model with lambda out performs model without
# This is the lambda model--- lambda on the slope and intercept---one cue (x)
# Testing different values lambda (lambda slope = lambda intercept): 0, 0.2,0.8

if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("lizzie", getwd())>0)) {
    setwd("~/Documents/git/teaching/stan/pmm/analyses")
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

### Flags to run models or just load data
runmodels <- FALSE

if(runmodels){
### Simulate data with and without lambda
nspecies = 40
nind = 10

nind <- rep(10, nspecies)
# Or... make it more uneven ... 
nind <- round(runif(nspecies, 1, 40))
# hist(nind, xlab="Simulated uneven sampling: n per sp", main="", breaks=100)

ndat <- sum(nind)

# Simulate species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Now set up the trait parameters
nruns <- 10
nparamL <- 8 # number of parameters including lp__
ntotL <- nparamL+2*nspecies

nparamN <- 6 # no lambda parameters
ntotN <- nparamN+2*nspecies
lambda = c(0, 0.2 , 1)
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

# Generate one data set
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsbf"]])         
slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)

dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind[i]),
                       intercept = rep(intercepts[i], nind[i]),
                       x1 = rnorm(n = nind[i], mean = 10, sd = 3),
                       trait1 = rep(slopes_bf[i], nind[i])
                       )
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
#sumDatLambOut <- rbind(sumDatLambOut, sumerLdf)


sumerLdf$true <- c(  sigma_y = param[["sigma_y"]],
                          lam_interceptsa = param[["lam_interceptsa"]],
                          sigma_interceptsa = param[["sigma_interceptsa"]],
                          lam_interceptsb = param[["lam_interceptsbf"]],
                          sigma_interceptsb = param[["sigma_interceptsbf"]],
                          slopes = slopes_bf,
                          b_z = param[["b_zf"]],
                          intercepts = intercepts,
                          a_z = param[["a_z"]],
                          lp = NA)
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


# add column of the true values to the dataframe
sumerNdf$true <- c( sigma_y = param[["sigma_y"]],
                    a_z = param[["a_z"]],
                    intercepts = intercepts,
                    sigma_a = param[["sigma_interceptsa"]],
                    b_z = param[["b_zf"]],
                    sigma_b = param[["sigma_interceptsbf"]],
                    slopes = slopes_bf,
                    lp = NA)

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
write.csv(sumDatL, paste("output/mdlOutLambdaRepped.csv", sep = ""))

sumDatN <- cbind(sumDatN, sumDatNoOut)
write.csv(sumDatN, paste("output/mdlOutNoLambdaRepped.csv", sep = ""))
# sumer <- data.frame(summary(test_new)$summary[c("a_z","lam_interceptsa","sigma_interceptsa", "b_z","lam_interceptsb" ,"sigma_interceptsb","sigma_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
# sumer <- data.frame(sumer)
# sumer$param <- do.call(rbind.data.frame, t(param))
# colnames(sumer)[colnames(sumer) == "c.4..0..0.2..0.6..0..0.1..0.01."] <- "testValue"
}

##########################################################################################
if(!runmodels){
# Pick if you want to run the even or uneven sampled sims
runeven <- FALSE
nspecies <- 40

if(runeven){ # highest lambda is 0.8
  sumDatL <- read.csv("output/mdlOutLambdaRepped.csv")
  sumDatN <- read.csv("output/mdlOutNoLambdaRepped.csv")
}

if(!runeven){ # highest lambda is 1
sumDatL <- read.csv("output/mdlOutLambdaUneven.csv")
sumDatN <- read.csv("output/mdlOutNoLambdaUneven.csv")
}


sumDatL$paramName <- c(  "sigma_y",
                         "lam_interceptsa" ,
                         "sigma_interceptsa",
                         "lam_interceptsb",
                         "sigma_interceptsb",
                         slopes = paste("b", seq(1:nspecies), sep = "_"),
                         "b_z",
                         intercepts = paste("a", seq(1:nspecies), sep = "_"),
                         "a_z"," lp__")


sumDatN$paramName <- c(  "sigma_y",
                         "mu_a" ,intercepts = paste("a", seq(1:nspecies), sep = "_"),
                         "sigma_a",
                         "mu_b",
                         slopes = paste("b", seq(1:nspecies), sep = "_"),
                         "sigma_b",
                         "lp__")
# Plotting model output:
  # 1. Compare the sp estimates from the test data to the sp estimates from the model, extract the coefficients and Rsq

# lambda model, pooled all values of lambda
slopesLamb <- sumDatL[grep("b\\[", sumDatL$X),]
slopesLamblam1 <- subset(slopesLamb, lambda==1|lambda==0.8)

mdlLamb <- summary(lm(true~mean, data = slopesLamb))

mdlLamb["r.squared"]$r.squared # 0.9998574

par(mfrow=c(2,2))
plot(slopesLamb$true ~ slopesLamb$mean, 
  xlab="estimated", ylab="true", main="lambda model: all data")
abline(lm(true~mean, data = slopesLamb))

plot(slopesLamblam1$true ~ slopesLamblam1$mean,
  xlab="estimated", ylab="true", main="lambda model: lamba=1 (or 0.8)")
abline(lm(true~mean, data = slopesLamblam1))


# NO lambda model
slopesNone <- sumDatN[grep("b\\[", sumDatN$X),]
slopesNonelam1 <- subset(slopesNone, lambda==1|lambda==0.8)

mdlNone <- summary(lm(true~mean, data = slopesNone))

mdlNone["r.squared"]$r.squared # 0.999858

plot(slopesNone$true ~ slopesNone$mean,
    xlab="estimated", ylab="true", main="hier. model: all data")
abline(lm(true~mean, data = slopesNone))

plot(slopesNonelam1$true ~ slopesNonelam1$mean,
    xlab="estimated", ylab="true", main="hier. model: lamba=1 (or 0.8)")
abline(lm(true~mean, data = slopesNonelam1))

  # 2. Calculate the diff of the mean value from the true value for the parameters

sumDatL$diffTrue <- sumDatL$mean - sumDatL$true


ggplot(sumDatL, aes(x= paramName, y=diffTrue)) +
  geom_boxplot(aes(col = as.factor(lambda)))

# 3. Count up when true value is in 50% interval
ifelse(sumDatL$X25.[1]<sumDatL$true[1] & sumDatL$X75.[1]>sumDatL$true[1], 1, 0)
sumDatL$truewithin50interval <- ifelse(sumDatL$X25.<sumDatL$true & sumDatL$X75.>sumDatL$true, 1, 0)
sumDatN$truewithin50interval <- ifelse(sumDatN$X25.<sumDatN$true & sumDatN$X75.>sumDatN$true, 1, 0)
sum(sumDatL$truewithin50interval, na.rm=TRUE)
sum(sumDatN$truewithin50interval, na.rm=TRUE)

# Do I have my signs right? Probably (otherwise these sums would go DOWN, not up)
sum(ifelse(sumDatL$X2.5.<sumDatL$true & sumDatL$X97.5.>sumDatL$true, 1, 0), na.rm=TRUE)
sum(ifelse(sumDatN$X2.5.<sumDatN$true & sumDatN$X97.5.>sumDatN$true, 1, 0), na.rm=TRUE)

# And now do it just for the species-level slopes ... 
sppDatL <- data.frame(lambda=character(), run=numeric(), intabstrue=numeric(),
  intntrue=numeric(), slopeabstrue=numeric(), slopentrue=numeric())
for(lamby in c(1:length(unique(sumDatL$lambda)))){
  subby <- subset(sumDatL, lambda==unique(sumDatL$lambda)[lamby])
  for(runhere in c(1:length(unique(subby$runs)))){
    runsubby <-  subset(subby, runs==unique(runs)[runhere])
    subbyint <- runsubby[grep("a\\[", runsubby$X),]
    subbyslope <- runsubby[grep("b\\[", runsubby$X),]
    dfhereadd <- data.frame(lambda=unique(sumDatL$lambda)[lamby], run=runhere, 
      intabstrue=sum(abs(subbyint$diffTrue), na.rm=TRUE), 
      intntrue=sum(subbyint$truewithin50interval, na.rm=TRUE), 
      slopeabstrue=sum(abs(subbyslope$diffTrue), na.rm=TRUE), 
      slopentrue=sum(subbyslope$truewithin50interval, na.rm=TRUE))
    sppDatL <- rbind(sppDatL, dfhereadd)
}
}

# And now do it just for the species-level slopes ... 
sppDatN <- data.frame(lambda=character(), run=numeric(), intabstrue=numeric(),
  intntrue=numeric(), slopeabstrue=numeric(), slopentrue=numeric())
for(lamby in c(1:length(unique(sumDatN$lambda)))){
  subby <- subset(sumDatN, lambda==unique(sumDatN$lambda)[lamby])
  for(runhere in c(1:length(unique(subby$runs)))){
    runsubby <-  subset(subby, runs==unique(runs)[runhere])
    subbyint <- runsubby[grep("a_sp\\[", runsubby$X),]
    subbyslope <- runsubby[grep("b\\[", runsubby$X),]
    dfhereadd <- data.frame(lambda=unique(sumDatN$lambda)[lamby], run=runhere, 
      intabstrue=sum(abs(subbyint$mean-subbyint$true), na.rm=TRUE), 
      intntrue=sum(subbyint$truewithin50interval, na.rm=TRUE), 
      slopeabstrue=sum(abs(subbyslope$mean-subbyslope$true), na.rm=TRUE), 
      slopentrue=sum(subbyslope$truewithin50interval, na.rm=TRUE))
    sppDatN <- rbind(sppDatN, dfhereadd)
  }
}

# Smaller is better (less deviations from true value)
aggregate(sppDatN["intabstrue"], sppDatN["lambda"], FUN=sum)
aggregate(sppDatL["intabstrue"], sppDatL["lambda"], FUN=sum)

aggregate(sppDatN["slopeabstrue"], sppDatN["lambda"], FUN=sum)
aggregate(sppDatL["slopeabstrue"], sppDatL["lambda"], FUN=sum)


# Bigger is better (50% interval captures true value more often)
aggregate(sppDatN["intntrue"], sppDatN["lambda"], FUN=sum)
aggregate(sppDatL["intntrue"], sppDatL["lambda"], FUN=sum)

aggregate(sppDatN["slopentrue"], sppDatN["lambda"], FUN=sum)
aggregate(sppDatL["slopentrue"], sppDatL["lambda"], FUN=sum)

  }
