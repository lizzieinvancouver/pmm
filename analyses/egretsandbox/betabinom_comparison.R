# Feb 2, 2025
# Started by V Van der Meersch
# Compare beta and binomial likelihoods!

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
require(progressr)

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
handlers(global = TRUE)
handlers("txtprogressbar")
set.seed(123456789)

source("functions/toolbox.R")
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")


run_models <- TRUE


# simulated data from prior distributions
nspecies <- 20 # number of species
spetree <- pbtree(n=nspecies, nsim=1, b=1, complete=FALSE,scale=1) # pure birth model
spetree$tip.label <- paste("s", 1:nspecies, sep="")
palette <- colorRampPalette(kippenberger)(nspecies)

nplots <- 121
simulated_data <- lapply(1:nplots, function(i){
  simulate_data_improved(spetree)
})

# fit models
if(run_models){
  sm <- stan_model("stan/binomiallikelihood.stan")
  plan(multisession, workers = 4) # here we use 4x4 = 16 cores (i.e. nested parallelization)
  datafit <- with_progress({
    p <- progressor(along = 1:nplots)
    future_lapply(simulated_data, function(simulated){
      
      mdl.data <- list(y = simulated$data$y,
                       N = nrow(simulated$data),
                       Nsp = nspecies,
                       sp = simulated$data$species,
                       x = simulated$data$x,
                       Ntrials = simulated$data$nseeds,
                       Vphy= simulated$vphy)
      
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
      
      sampler_params  <- get_sampler_params(fit, inc_warmup = FALSE)
      diagnostics <- list(
        max_treedepth= max(sapply(sampler_params, function(x) max(x[, "treedepth__"]))),
        max_divergence = max(sapply(sampler_params, function(x) sum(x[, "divergent__"]))),
        max_rhat = max(summ$Rhat, na.rm = TRUE),
        min_ess = min(summ$n_eff, na.rm = TRUE)
      )
      
      estimates <- get_estimates(summ)
      
      posteriors <- c(extract(fit, pars = c("lambda_a", "lambda_b")), chain = list(rep(1:4, each = 1000)), limits = list(c(0,1)))
      
      p()
      
      return(list(sim = simulated$data, 
                  est = estimates,
                  post = posteriors,
                  par = simulated$params,
                  diag = diagnostics))
    }, future.seed=TRUE)})
  plan(sequential);gc()
  saveRDS(datafit, file = "output/comparison/bionomiallikelihood.rds")
}else{
  datafit <- readRDS(file = "output/comparison/bionomiallikelihood.rds")
}

# order next plots by sigma_b values
sigmasb <- sapply(1:nplots, function(i) datafit[[i]][["par"]]$sigma_b)

# plot simulated data (binomial mode)
pdf(file=paste0("figures/comparison/simulated_data_count.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  diagnostics <- datafit[[i]][["diag"]]
  warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 100)
  
  plot.new()
  plot.window(xlim = c(-3,3), ylim = c(0,100))
  if(warning){
    rect(par("usr")[1], par("usr")[3],
         par("usr")[2], par("usr")[4],
         col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  }
  grid()
  points(x = simulated_data[[i]]$data$x, y =  simulated_data[[i]]$data$y, 
         pch = 19, cex = 0.2,
         col = palette[simulated_data[[i]]$data$species])
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", simulated_data[[i]]$params$a_z, ", l_a=", simulated_data[[i]]$params$lambda_a, ", s_a=", simulated_data[[i]]$params$sigma_a, "\n",
               "b_z=", simulated_data$params$b_z, ", l_b=", simulated_data[[i]]$params$lambda_b, ", s_b=", simulated_data[[i]]$params$sigma_b), 
        adj=0, cex.main = 0.5)
  
  
}
dev.off()

# plot simulated data (beta model)
pdf(file=paste0("figures/comparison/simulated_data_percent.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  diagnostics <- datafit[[i]][["diag"]]
  warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 100)
  
  plot.new()
  plot.window(xlim = c(-3,3), ylim = c(0,1))
  if(warning){
    rect(par("usr")[1], par("usr")[3],
         par("usr")[2], par("usr")[4],
         col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  }
  grid()
  points(x = simulated_data[[i]]$data$x, y =  simulated_data[[i]]$data$y/simulated_data[[i]]$data$nseeds, 
         pch = 19, cex = 0.2,
         col = palette[simulated_data[[i]]$data$species])
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", simulated_data[[i]]$params$a_z, ", l_a=", simulated_data[[i]]$params$lambda_a, ", s_a=", simulated_data[[i]]$params$sigma_a, "\n",
               "b_z=", simulated_data$params$b_z, ", l_b=", simulated_data[[i]]$params$lambda_b, ", s_b=", simulated_data[[i]]$params$sigma_b), 
        adj=0, cex.main = 0.5)
  
  
}
dev.off()

# plot estimated intercepts
pdf(file=paste0("figures/comparison/intercepts_binomial.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  estimates <- datafit[[i]][["est"]]
  params <- datafit[[i]][["par"]]
  diagnostics <- datafit[[i]][["diag"]]
  
  warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 100)
  
  plot.new()
  limits <- c(-5,5)
  plot.window(xlim = limits, ylim = limits)
  if(warning){
    rect(par("usr")[1], par("usr")[3],
         par("usr")[2], par("usr")[4],
         col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  }
  grid()
  abline(a=0, b=1, col = "grey80")
  segments(y0 = estimates$intercepts.q2.5, x0 = params$intercepts, 
           y1 = estimates$intercepts.q97.5, x1 = params$intercepts,
           col = palette)
  points(x = params$intercepts, y =  estimates$intercepts.mean, 
         pch = 19, cex = 0.5, col = palette)
  title(ylab = "Estimated", xlab = "Observed", cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
               "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b), 
        adj=0, cex.main = 0.5)
  text(x = -4.7, y = 4, label = paste("max_tdeepth=", diagnostics$max_treedepth), cex = 0.45, 
       col = ifelse(diagnostics$max_treedepth < 10, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 3.5, label = paste("max_dvg=", diagnostics$max_divergence), cex = 0.45, 
       col = ifelse(diagnostics$max_divergence == 0, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 3, label = paste("max_rhat=", round(diagnostics$max_rhat,2)), cex = 0.45, 
       col = ifelse(diagnostics$max_rhat < 1.04, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 2.5, label = paste("min_ess=", round(diagnostics$min_ess),1), cex = 0.45, 
       col = ifelse(diagnostics$min_ess > 100, "grey40", "darkred"), adj=0) 
  
}
dev.off()



