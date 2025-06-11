
rm(list = ls())
library(ape)
library(geiger)
library(phytools)
library(rstan)
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

nspecies <- 40 # let's say we have 40 species 

# stochastic birth-death trees (from Deirdre)
maxlength <- runif(1,10,20)
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=maxlength, type = "continuous") 
# with a total length bewteen 10 and 20 (scale = maxlength)
# note: we could set a length = 1, but we want to avoid any confusion 
# between the covariance matrix [vcv(spetree, corr = FALSE)]
# and the correlation matrix [vcv(spetree, corr = TRUE)]...

spetree$tip.label <- paste("s", 1:nspecies, sep="")
palette <- colorRampPalette(kippenberger)(nspecies)


plot(spetree, show.tip.label = FALSE, y.lim = c(0,45), x.lim = c(0,maxlength))
segments(x0 = 0, y0 = 43, x1 = maxlength, col = "#AE2565FF", lwd=1.5)
text(x = maxlength/2, y = 44.2, label = expression(paste("Tree length = ", varphi)), cex = 1,  col = "#AE2565FF")

# distance between species
distance <- cophenetic(spetree)

# phylo covariance matrix
covariance_matrix <- vcv(spetree, corr = FALSE)
correlation_matrix <- vcv(spetree, corr = TRUE)

dataplot <- unique(data.frame(
  dist = c(distance),
  cov = c(covariance_matrix),
  corr = c(correlation_matrix)
))

par(mfrow = c(1,2))

plot.new()
plot.window(xlim = c(0,max(distance)), ylim = c(0,maxlength))
grid()
axis(1, las = 1, cex.axis = 1, tck=-0.02, at = c(0,max(distance)), labels=c("0", expression(paste("2", varphi))))
axis(2, las = 1, cex.axis = 1, tck=-0.02, at = c(0,maxlength), labels=c("0", expression(paste(varphi))))
title(main = "Covariogram", xlab = expression(paste(Delta,x)),  line=1, cex.main = 1)
lines(dataplot$dist, dataplot$cov, col="white", lwd=4)
lines(dataplot$dist, dataplot$cov, col="#AE2565FF", lwd=2)

plot.new()
plot.window(xlim = c(0,max(distance)), ylim = c(0,1))
grid()
axis(1, las = 1, cex.axis = 1, tck=-0.02, at = c(0,max(distance)), labels=c("0", expression(paste("2", varphi))))
axis(2, las = 1, cex.axis = 1, tck=-0.02, at = c(0,1), labels=c("0", "1"))
title(main = "Correlogram", xlab = expression(paste(Delta,x)),  line=1, cex.main = 1)
lines(dataplot$dist, dataplot$corr, col="white", lwd=4)
lines(dataplot$dist, dataplot$corr, col="#AE2565FF", lwd=2)


# parameters
params <- list(a_z = rnorm(1,0,1.5), # root value intercept
               lambda_a = rbeta(1,1.5,1.5), # lambda intercept
               sigma_a = abs(rnorm(1,0,1)) # rate of evolution intercept
)

params$lambda_a <- 0.7

## Simulate some data ## 
# phylogenetic structure (from Deirdre)
scaledtree_intercept <- rescale(spetree, model = "lambda", params$lambda_a)         
params$intercepts <- fastBM(scaledtree_intercept,
                            a = params$a_z, mu = 0, 
                            sig2 = params$sigma_a ^ 2)

distance <- cophenetic(spetree)
covariance_matrix <- vcv(scaledtree_intercept, corr = FALSE)
dataplot <- unique(data.frame(
  dist = c(distance),
  cov = c(covariance_matrix)
))
dataplot <- dataplot[order(dataplot$dist),]

dataplot_sub <- dataplot[dataplot$dist != 0,]
plot.new()
plot.window(xlim = c(0,max(distance)), ylim = c(0,maxlength))
grid()
axis(1, las = 1, cex.axis = 1, tck=-0.02, at = c(0,max(distance)), labels=c("0", expression(paste("2", varphi))))
axis(2, las = 1, cex.axis = 1, tck=-0.02, at = c(0,params$lambda_a*maxlength,maxlength), labels=c("0", expression(paste(lambda, "*", varphi)), expression(paste(varphi))))
title(main = "Covariogram", xlab = expression(paste(Delta,x)),  line=1, cex.main = 1)
lines(dataplot_sub$dist, dataplot_sub$cov, col="white", lwd=4)
lines(dataplot_sub$dist, dataplot_sub$cov, col="#AE2565FF", lwd=2)
points(dataplot[dataplot$dist == 0, "dist"], dataplot[dataplot$dist == 0, "cov"], 
       col="#AE2565FF", pch = 20)



# We consider that a study looks only  at one species
# and that a species might be considered by several studies
nstudies_perspecies <- round(runif(nspecies, 1, 3)) # no. of different studies per species
nstudies <- sum(nstudies_perspecies) # total number of studies
species_study <- rep(1:nspecies, times = nstudies_perspecies) # species id considered by each study
ntreat_perstudies <- round(runif(nstudies, 2, 6)) # no. of different treatments applied per study
nexps <- sum(ntreat_perstudies) # total no. of experiments (i.e. studies*treatments)
ntrialseeds_perexp <- round(runif(nexps, 10, 100)) # no. of seeds per experiments
studies <- rep(1:nstudies, times = ntreat_perstudies)
species <- rep(species_study, times = ntreat_perstudies) # species id per experiments

# experimental observations
# x <- rnorm(n = nexps, mean = 0, sd = 1) # treatments applied
yhat <- plogis(params$intercepts[species], 
               location = 0, scale = 1)
y = rbinom(n = nexps, size = ntrialseeds_perexp, prob = yhat)

simulated <- list(
  params = lapply(params,round,2), 
  data = data.frame(species, studies, nseeds = ntrialseeds_perexp, y),
  vphy = vcv(spetree, corr = TRUE))

smbin <- stan_model("/home/victor/projects/phylo/binomiallikelihood.stan")

mdl.data <- list(y = simulated$data$y,
                 N = nrow(simulated$data),
                 Nsp = nspecies,
                 sp = simulated$data$species,
                 # x = simulated$data$x,
                 Ntrials = simulated$data$nseeds,
                 Vphy= simulated$vphy)

fit <- sampling(smbin, mdl.data, 
                iter = 4000, warmup = 3000,
                chains = 4, cores = 4)

