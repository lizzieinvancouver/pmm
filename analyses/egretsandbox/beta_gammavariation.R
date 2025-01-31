# Jan 27, 2025
# Started by V Van der Meersch
# Largely inspired by Deirdre work



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

## Simulate some data ##
nspecies <- 40
nind_species <- round(runif(nspecies, 1, 50))
n <- sum(nind_species)

# species tree with a pure birth model
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nspecies, sep="")

params <- list(a_z = -1, # root value intercept
              lambda_a = 0.8, # lambda intercept
              sigma_a = 0.5, # rate of evolution intercept
              b_z = 0, # root value trait1 slope
              lambda_b = 0.1, # lambda trait1
              sigma_b = 0.2, # rate of evolution trait1
              gamma = NA # intercept dispersion parameter
)

## Vary the gamma (~ dispersion parameter)
gammas <- c(10,5,1,0.05)
simulated_data <- lapply(gammas, function(g){
  params$gamma <- g
  list(data = simulate_data(spetree, nind_species, params), spetree = spetree)
})

sm <- stan_model("beta.stan")


plan(multisession, workers = 4) # here we use 4x4 = 16 cores (i.e. nested parallelization)
datafit <- future_lapply(simulated_data, function(simulated){
  
  mdl.data <- list(y = simulated$data$y,
                   N = nrow(simulated$data),
                   Nsp = nspecies,
                   sp = simulated$data$sp,
                   x = simulated$data$x,
                   Vphy= vcv(simulated$spetree, corr = TRUE))
  
  # fit <- stan("beta.stan",
  #             data = mdl.data,
  #             iter = 4000,
  #             warmup = 3000,
  #             chains = 4,
  #             cores = 4)
  
  fit <- sampling(sm, mdl.data, 
                  iter = 4000, warmup = 3000,
                  chains = 4, cores = 4)
  
  summ <- data.frame(summary(fit)[["summary"]])
  
  estimates <- data.frame(sigma_y.mean = summ["gamma" , "mean"],
                           sigma_y.se = summ["gamma" , "se_mean"],
                           lambda_a.mean = summ["lambda_a" , "mean"],
                           lambda_a.se = summ["lambda_a" , "se_mean"],
                           sigma_a.mean = summ["sigma_a" , "mean"],
                           sigma_a.se = summ["sigma_a" , "se_mean"],
                           lambda_b.mean = summ["lambda_b" , "mean"],
                           lambda_b.se = summ["lambda_b" , "se_mean"],
                           sigma_b.mean = summ["sigma_b" , "mean"],
                           sigma_b.se = summ["sigma_b" , "se_mean"],
                           slopes.mean = summ[paste0("b[", 1:nspecies, "]"), "mean"],
                           slopes.q2.5 = summ[paste0("b[", 1:nspecies, "]"), "X2.5."],
                           slopes.q97.5 = summ[paste0("b[", 1:nspecies, "]"), "X97.5."],
                           intercepts.mean = summ[paste0("a[", 1:nspecies, "]"), "mean"],
                           intercepts.q2.5 = summ[paste0("a[", 1:nspecies, "]"), "X2.5."],
                           intercepts.q97.5 = summ[paste0("a[", 1:nspecies, "]"), "X97.5."],
                           a_z.mean = summ["a_z", "mean"],
                           a_z.q2.5 = summ["a_z", "X2.5."],
                           a_z.q97.5 = summ["a_z", "X97.5."])
  
  posteriors <- c(extract(fit, pars = c("lambda_a", "lambda_b")), chain = list(rep(1:4, each = 1000)), limits = list(c(0,1)))
  
  observed <- data.frame(unique(simulated$data[c("sp", "slopes", "intercepts")]), par1 = params$lambda_a, par2 = params$lambda_b)
  
  return(list(sim = simulated$data, 
              obs = observed, 
              est = estimates,
              post = posteriors,
              par = params))
}, future.seed=TRUE)
plan(sequential);gc()

pdf(file=paste0("figures/gamma_variation_broot",params$b_z,".pdf"))
par(mfrow = c(length(gammas),4),
    mex=0.6, mar=c(5,4,3,2)+.1)
for(i in 1:length(gammas)){
  plot_stuff(datafit[[i]][["est"]], datafit[[i]][["obs"]], datafit[[i]][["sim"]], datafit[[i]][["post"]], comments = paste0(', disp.=', gammas[i]))
}
dev.off()

