# Jan 29, 2025
# Started by V Van der Meersch

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pmm/analyses")
} else if(length(grep("lizzie", getwd())>0)) {
  setwd("~/Documents/git/teaching/stan/pmm/analyses")
} else if(length(grep("victor", getwd())>0)) {
  setwd("/home/victor/projects/pmm/analyses/egretsandbox")
} else{
  setwd("/home/deirdre/pmm") # for midge
}

# Load packages
require(ape)
require(geiger)
require(phytools)
require(rstan)
require(shinystan)
require(future.apply)

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(1234567)

source("functions/toolbox.R")
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")


## Simulate some data ##
nspecies <- 40
nind_species <- round(runif(nspecies, 1, 50))
n <- sum(nind_species)
palette <- colorRampPalette(kippenberger)(nspecies)

# species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

# Priors in beta.stan
# a_z ~ normal(0, 1.5); 
# b_z ~ normal(0.5, 1); 
# 
# lambda_a ~ beta(1, 1);
# sigma_a ~ normal(0, 1);
# 
# lambda_b ~ beta(1, 1);
# sigma_b ~ normal(0, 1);
# 
# gamma ~ cauchy(0, 3); 

nplots <- 121

pdf(file=paste0("figures/priors_check2.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in 1:nplots){
  
  params <- list(a_z = rnorm(1,0,1.5), # root value intercept
                 lambda_a = rbeta(1,1,1), # lambda intercept
                 sigma_a = rnorm(1,0,1), # rate of evolution intercept
                 b_z = abs(rnorm(1,0.5,1)), # root value trait1 slope
                 lambda_b = rbeta(1,1,1), # lambda trait1
                 sigma_b = abs(rnorm(1,0,1)), # rate of evolution trait1
                 gamma = abs(rcauchy(1,0,3)) # intercept dispersion parameter
  )
  
  simulated_data <- simulate_data(spetree, nind_species, params,
                                  mux = 0, sigmax = 2)
  
  params <- lapply(params, round, 2)
  
  plot.new()
  plot.window(xlim = c(-5,5), ylim = c(0,1))
  grid()
  points(x = simulated_data$x, y =  simulated_data$y, 
         pch = 19, cex = 0.2,
         col = palette[simulated_data$sp])
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
                      "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b, ", g=", params$gamma), 
        adj=0, cex.main = 0.5)
  
  
}
dev.off()

