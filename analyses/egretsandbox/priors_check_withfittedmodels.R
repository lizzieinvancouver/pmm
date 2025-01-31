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
require(progressr)

rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
handlers(global = TRUE)
handlers("txtprogressbar")
set.seed(1234567)

source("functions/toolbox.R")
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

run_models <- TRUE


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
# lambda_a ~ beta(1.5, 1.5);
# sigma_a ~ normal(0, 1);
# 
# lambda_b ~ beta(1.5, 1.5);
# sigma_b ~ normal(0, 1);
# 
# // phi ~ cauchy(0, 3); previously
# phi ~ gamma(4, 0.1) // give more plausible shapes

nplots <- 121

# simulated data from prior distributions
simulated_data <- lapply(1:nplots, function(i){
  
  params <- list(a_z = rnorm(1,0,1.5), # root value intercept
                 lambda_a = rbeta(1,1.5,1.5), # lambda intercept
                 sigma_a = abs(rnorm(1,0,1)), # rate of evolution intercept
                 b_z = rnorm(1,0.5,1), # root value trait1 slope
                 lambda_b = rbeta(1,1.5,1.5), # lambda trait1
                 sigma_b = abs(rnorm(1,0,1)), # rate of evolution trait1
                 # phi = abs(rcauchy(1,0,3)) # dispersion parameter
                 phi = rgamma(1,4,0.1)
  )
  params <- lapply(params, round, 2)
  list(data = simulate_data(spetree, nind_species, params), spetree = spetree, params = params)
  
})

# fit models
if(run_models){
  sm <- stan_model("stan/betalikelihood.stan")
  plan(multisession, workers = 4) # here we use 4x4 = 16 cores (i.e. nested parallelization)
  datafit <- with_progress({
    p <- progressor(along = 1:nplots)
    future_lapply(simulated_data, function(simulated){
    
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
    
    sampler_params  <- get_sampler_params(fit, inc_warmup = FALSE)
    diagnostics <- list(
      max_treedepth= max(sapply(sampler_params, function(x) max(x[, "treedepth__"]))),
      max_divergence = max(sapply(sampler_params, function(x) sum(x[, "divergent__"]))),
      max_rhat = max(summ$Rhat, na.rm = TRUE),
      min_ess = min(summ$n_eff, na.rm = TRUE)
    )
    
    estimates <- data.frame(phi.mean = summ["phi" , "mean"],
                            phi.q2.5 = summ["phi" , "X2.5."],
                            phi.q97.5 = summ["phi" , "X97.5."],
                            lambda_a.mean = summ["lambda_a" , "mean"],
                            lambda_a.q2.5 = summ["lambda_a" , "X2.5."],
                            lambda_a.q97.5  = summ["lambda_a" , "X97.5."],
                            sigma_a.mean = summ["sigma_a" , "mean"],
                            sigma_a.q2.5 = summ["sigma_a" , "X2.5."],
                            sigma_a.q97.5  = summ["sigma_a" , "X97.5."],
                            lambda_b.mean = summ["lambda_b" , "mean"],
                            lambda_b.q2.5 = summ["lambda_b" , "X2.5."],
                            lambda_b.q97.5  = summ["lambda_b" , "X97.5."],
                            sigma_b.mean = summ["sigma_b" , "mean"],
                            sigma_b.q2.5 = summ["sigma_b" , "X2.5."],
                            sigma_b.q97.5  = summ["sigma_b" , "X97.5."],
                            slopes.mean = summ[paste0("b[", 1:nspecies, "]"), "mean"],
                            slopes.q2.5 = summ[paste0("b[", 1:nspecies, "]"), "X2.5."],
                            slopes.q97.5 = summ[paste0("b[", 1:nspecies, "]"), "X97.5."],
                            intercepts.mean = summ[paste0("a[", 1:nspecies, "]"), "mean"],
                            intercepts.q2.5 = summ[paste0("a[", 1:nspecies, "]"), "X2.5."],
                            intercepts.q97.5 = summ[paste0("a[", 1:nspecies, "]"), "X97.5."],
                            a_z.mean = summ["a_z", "mean"],
                            a_z.q2.5 = summ["a_z", "X2.5."],
                            a_z.q97.5 = summ["a_z", "X97.5."],
                            b_z.mean = summ["b_z", "mean"],
                            b_z.q2.5 = summ["b_z", "X2.5."],
                            b_z.q97.5 = summ["b_z", "X97.5."])
    
    posteriors <- c(extract(fit, pars = c("lambda_a", "lambda_b")), chain = list(rep(1:4, each = 1000)), limits = list(c(0,1)))
    
    observed <- data.frame(unique(simulated$data[c("sp", "slopes", "intercepts")]), par1 = simulated$params$lambda_a, par2 = simulated$params$lambda_b)
    
    p()
    
    return(list(sim = simulated$data, 
                obs = observed, 
                est = estimates,
                post = posteriors,
                par = simulated$params,
                diag = diagnostics))
  }, future.seed=TRUE)})
  plan(sequential);gc()
  saveRDS(datafit, file = "output/models_simulateddata_frompriors.rds")
}else{
  datafit <- readRDS(file = "output/models_simulateddata_frompriors.rds")
}

# order next plots by phi values
# phis <- sapply(1:nplots, function(i) datafit[[i]][["par"]]$phi)
# order next plots by sigma_b values
sigmasb <- sapply(1:nplots, function(i) datafit[[i]][["par"]]$sigma_b)

# plot simulated data
pdf(file=paste0("figures/priors_check.pdf"), height = 15, width = 18)
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
  points(x = simulated_data[[i]]$data$x, y =  simulated_data[[i]]$data$y, 
         pch = 19, cex = 0.2,
         col = palette[simulated_data[[i]]$data$sp])
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", simulated_data[[i]]$params$a_z, ", l_a=", simulated_data[[i]]$params$lambda_a, ", s_a=", simulated_data[[i]]$params$sigma_a, "\n",
               "b_z=", simulated_data[[i]]$params$b_z, ", l_b=", simulated_data[[i]]$params$lambda_b, ", s_b=", simulated_data[[i]]$params$sigma_b, ", g=", simulated_data[[i]]$params$phi), 
        adj=0, cex.main = 0.5)
  
  
}
dev.off()

# plot estimated intercepts
pdf(file=paste0("figures/priors_check_wfitintercept.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  estimates <- datafit[[i]][["est"]]
  observed <- datafit[[i]][["obs"]]
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
  segments(y0 = estimates$intercepts.q2.5, x0 = observed$intercepts, 
           y1 = estimates$intercepts.q97.5, x1 = observed$intercepts,
           col = palette)
  points(x = observed$intercepts, y =  estimates$intercepts.mean, 
         pch = 19, cex = 0.5, col = palette)
  title(ylab = "Estimated", xlab = "Observed", cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
               "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b, ", g=", params$phi), 
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

# plot estimated slopes
pdf(file=paste0("figures/priors_check_wfitslopes.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  estimates <- datafit[[i]][["est"]]
  observed <- datafit[[i]][["obs"]]
  params <- datafit[[i]][["par"]]
  diagnostics <- datafit[[i]][["diag"]]
  
  warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 400)
  
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
  segments(y0 = estimates$slopes.q2.5, x0 = observed$slopes, 
           y1 = estimates$slopes.q97.5, x1 = observed$slopes,
           col = palette)
  points(x = observed$slopes, y =  estimates$slopes.mean, 
         pch = 19, cex = 0.5, col = palette)
  title(ylab = "Estimated", xlab = "Observed", cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
               "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b, ", g=", params$phi), 
        adj=0, cex.main = 0.5)
  text(x = -4.7, y = 4, label = paste("max_tdeepth=", diagnostics$max_treedepth), cex = 0.45, 
       col = ifelse(diagnostics$max_treedepth < 10, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 3.5, label = paste("max_dvg=", diagnostics$max_divergence), cex = 0.45, 
       col = ifelse(diagnostics$max_divergence == 0, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 3, label = paste("max_rhat=", round(diagnostics$max_rhat,2)), cex = 0.45, 
       col = ifelse(diagnostics$max_rhat < 1.04, "grey40", "darkred"), adj=0) 
  text(x = -4.7, y = 2.5, label = paste("min_ess=", round(diagnostics$min_ess),1), cex = 0.45, 
       col = ifelse(diagnostics$min_ess > 400, "grey40", "darkred"), adj=0) 
  
}
dev.off()


pdf(file=paste0("figures/priors_check_wglobalparams.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1.4,0,0,0)+0.7, mgp=c(0,0.5,0))
for(i in order(sigmasb)){
  
  estimates <- datafit[[i]][["est"]]
  observed <- datafit[[i]][["obs"]]
  params <- datafit[[i]][["par"]]
  diagnostics <- datafit[[i]][["diag"]]
  
  warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 400)
  
  plot.new()
  limits <- c(1-0.25,7+0.25)
  plot.window(xlim = limits, ylim = c(-5,5))
  if(warning){
    rect(par("usr")[1], par("usr")[3],
         par("usr")[2], par("usr")[4],
         col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  }
  grid()
  
  # phi, mean and std of prior (hack here...)
  mphi <- mean(rgamma(10000,4,0.1))
  stdphi <- sd(rgamma(10000,4,0.1))
  
  # root intercept
  est_df <- 
    data.frame(mean = c(estimates$a_z.mean[1], estimates$lambda_a.mean[1], estimates$sigma_a.mean[1], 
                        estimates$b_z.mean[1], estimates$lambda_b.mean[1], estimates$sigma_b.mean[1], (estimates$phi.mean[1]-mphi)/stdphi),
               q2.5 = c(estimates$a_z.q2.5[1], estimates$lambda_a.q2.5[1], estimates$sigma_a.q2.5[1], 
                        estimates$b_z.q2.5[1], estimates$lambda_b.q2.5[1], estimates$sigma_b.q2.5[1], (estimates$phi.q2.5[1]-mphi)/stdphi), 
               q97.5 = c(estimates$a_z.q97.5[1], estimates$lambda_a.q97.5[1], estimates$sigma_a.q97.5[1],
                         estimates$b_z.q97.5[1], estimates$lambda_b.q97.5[1], estimates$sigma_b.q97.5[1], (estimates$phi.q97.5[1]-mphi)/stdphi),
               x = c(1:7))
  labs <- c(expression(a [root]), expression(lambda  [a]), expression(sigma  [a]),
            expression(b [root]), expression(lambda  [b]), expression(sigma  [b]), expression(phi [scaled]))
  
  param_df <- 
    data.frame(mean = c(params$a_z, params$lambda_a, params$sigma_a, params$b_z, params$lambda_b, params$sigma_b, (params$phi-mphi)/stdphi), 
               x= c(1:7))
  
  segments(y0 = param_df$mean, x0 = param_df$x-0.25, 
           y1 = param_df$mean, x1 = param_df$x+0.25,
           col = "grey50")
  segments(y0 = est_df$q2.5, x0 = est_df$x, 
           y1 = est_df$q97.5, x1 = est_df$x,
           col = palette)
  points(x = est_df$x, y =  est_df$mean, 
         pch = 19, cex = 0.5, col = palette)
  axis(1, las = 2, cex.axis = 0.7, tck=-0.02, at = 1:7, labels = labs)
  axis(2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  # title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
  #              "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b, ", g=", params$gamma), 
  #       adj=0, cex.main = 0.5)
  text(x = 0.7, y = 4.5, label = paste("max_tdeepth=", diagnostics$max_treedepth), cex = 0.45, 
       col = ifelse(diagnostics$max_treedepth < 10, "grey40", "darkred"), adj=0) 
  text(x = 0.7, y = 4, label = paste("max_dvg=", diagnostics$max_divergence), cex = 0.45, 
       col = ifelse(diagnostics$max_divergence == 0, "grey40", "darkred"), adj=0) 
  text(x = 0.7, y = 3.5, label = paste("max_rhat=", round(diagnostics$max_rhat,2)), cex = 0.45, 
       col = ifelse(diagnostics$max_rhat < 1.04, "grey40", "darkred"), adj=0) 
  text(x = 0.7, y = 3, label = paste("min_ess=", round(diagnostics$min_ess),1), cex = 0.45, 
       col = ifelse(diagnostics$min_ess > 400, "grey40", "darkred"), adj=0) 
  
}
dev.off()



