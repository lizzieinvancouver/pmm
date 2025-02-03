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
  
  smbin <- stan_model("stan/binomiallikelihood.stan")
  plan(multisession, workers = 4) # here we use 4x4 = 16 cores (i.e. nested parallelization)
  databin <- with_progress({
    p <- progressor(along = 1:nplots)
    future_lapply(simulated_data, function(simulated){
      
      mdl.data <- list(y = simulated$data$y,
                       N = nrow(simulated$data),
                       Nsp = nspecies,
                       sp = simulated$data$species,
                       x = simulated$data$x,
                       Ntrials = simulated$data$nseeds,
                       Vphy= simulated$vphy)
    
      fit <- sampling(smbin, mdl.data, 
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
  saveRDS(databin, file = "output/comparison/bionomiallikelihood.rds")
  
  smordbeta <- stan_model("stan/orderedbetalikelihood.stan")
  plan(multisession, workers = 4) 
  databeta <- with_progress({
    p <- progressor(along = 1:nplots)
    future_lapply(simulated_data, function(simulated){
      
      simulated$data$yperc <- simulated$data$y/simulated$data$nseeds
      
      mdl.data <- list(N_degen=sum(simulated$data$yperc %in% c(0,1)),
                       N_prop=sum(simulated$data$yperc>0 & simulated$data$yperc<1),
                       
                       Nsp = nspecies,
                       sp_degen = array(simulated$data$species[simulated$data$yperc %in% c(0,1)],
                                        dim = sum(simulated$data$yperc %in% c(0,1))),
                       sp_prop = simulated$data$species[simulated$data$yperc>0 & simulated$data$yperc<1],
                       
                       
                       y_degen=array(simulated$data$yperc[simulated$data$yperc %in% c(0,1)],
                                     dim = sum(simulated$data$yperc %in% c(0,1))),
                       y_prop=simulated$data$yperc[simulated$data$yperc>0 & simulated$data$yperc<1],
                       
                       x_degen=array(simulated$data$x[simulated$data$yperc %in% c(0,1)],
                                     dim = sum(simulated$data$yperc %in% c(0,1))),
                       x_prop=simulated$data$x[simulated$data$yperc>0 & simulated$data$yperc<1],
                       
                       Vphy= simulated$vphy)
      
      fit <- sampling(smordbeta, mdl.data, 
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
  saveRDS(databeta, file = "output/comparison/betalikelihood.rds")
}else{
  databin <- readRDS(file = "output/comparison/bionomiallikelihood.rds")
  databeta <- readRDS(file = "output/comparison/betalikelihood.rds")
}

# order next plots by sigma_b values
sigmasb <- sapply(1:nplots, function(i) databin[[i]][["par"]]$sigma_b)

# plot simulated data (binomial mode)
pdf(file=paste0("figures/comparison/simulated_data_count.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  diagnostics <- databin[[i]][["diag"]]
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
pdf(file=paste0("figures/comparison/simulated_data_prop.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(i in order(sigmasb)){
  
  diagnostics <- databeta[[i]][["diag"]]
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

# compare estimated intercepts - by both models
pdf(file=paste0("figures/comparison/intercepts_comparison.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(0.7,0.4,0,0)+0.7, mgp=c(0.3,0,0))
for(i in order(sigmasb)){
  
  estimatesbin <- databin[[i]][["est"]]
  estimatesbeta <- databeta[[i]][["est"]]
  params <- databeta[[i]][["params"]]
  
  # warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 100)
  
  plot.new()
  limits <- c(-5,5)
  plot.window(xlim = limits, ylim = limits)
  # if(warning){
  #   rect(par("usr")[1], par("usr")[3],
  #        par("usr")[2], par("usr")[4],
  #        col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  # }
  grid()
  abline(a=0, b=1, col = "grey80")
  segments(y0 = estimatesbin$intercepts.q2.5, x0 = estimatesbeta$intercepts.mean, 
           y1 = estimatesbin$intercepts.q97.5, x1 = estimatesbeta$intercepts.mean,
           col = palette)
  segments(x0 = estimatesbeta$intercepts.q2.5, y0 = estimatesbin$intercepts.mean, 
           x1 = estimatesbeta$intercepts.q97.5, y1 = estimatesbin$intercepts.mean,
           col = palette)
  points(x = estimatesbeta$intercepts.mean, y =  estimatesbin$intercepts.mean, 
         pch = 19, cex = 0.5, col = palette)
  title(ylab = "Estimated (binomial model)", cex.lab = 0.5)
  title(xlab = "Estimated (beta model)", cex.lab = 0.5, line = -0.1)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
               "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b), 
        adj=0, cex.main = 0.5)
  
}
dev.off()

# compare estimated slopes - by both models
pdf(file=paste0("figures/comparison/slopes_comparison.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(0.7,0.4,0,0)+0.7, mgp=c(0.3,0,0))
for(i in order(sigmasb)){
  
  estimatesbin <- databin[[i]][["est"]]
  estimatesbeta <- databeta[[i]][["est"]]
  params <- databeta[[i]][["params"]]
  
  # warning <- (diagnostics$max_treedepth >= 10 | diagnostics$max_divergence > 0 | diagnostics$max_rhat >= 1.04 | diagnostics$min_ess <= 100)
  
  plot.new()
  limits <- c(-5,5)
  plot.window(xlim = limits, ylim = limits)
  # if(warning){
  #   rect(par("usr")[1], par("usr")[3],
  #        par("usr")[2], par("usr")[4],
  #        col = rgb(0.545,0,0, alpha = 0.1), border = NA)
  # }
  grid()
  abline(a=0, b=1, col = "grey80")
  segments(y0 = estimatesbin$slopes.q2.5, x0 = estimatesbeta$slopes.mean, 
           y1 = estimatesbin$slopes.q97.5, x1 = estimatesbeta$slopes.mean,
           col = palette)
  segments(x0 = estimatesbeta$slopes.q2.5, y0 = estimatesbin$slopes.mean, 
           x1 = estimatesbeta$slopes.q97.5, y1 = estimatesbin$slopes.mean,
           col = palette)
  points(x = estimatesbeta$slopes.mean, y =  estimatesbin$slopes.mean, 
         pch = 19, cex = 0.5, col = palette)
  title(ylab = "Estimated (binomial model)", cex.lab = 0.5)
  title(xlab = "Estimated (beta model)", cex.lab = 0.5, line = -0.1)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  title(paste0("a_z=", params$a_z, ", l_a=", params$lambda_a, ", s_a=", params$sigma_a, "\n",
               "b_z=", params$b_z, ", l_b=", params$lambda_b, ", s_b=", params$sigma_b), 
        adj=0, cex.main = 0.5)
  
}
dev.off()

# Compute mean/std of priors (to scale estimates in the following)
priors <- list(a_z = rnorm(1e6,0,1.5), # root value intercept
               lambda_a = rbeta(1e6,1.5,1.5), # lambda intercept
               sigma_a = abs(rnorm(1e6,0,1)), # rate of evolution intercept
               b_z = rnorm(1e6,0.5,1), # root value trait1 slope
               lambda_b = rbeta(1e6,1.5,1.5), # lambda trait1
               sigma_b = abs(rnorm(1e6,0,1)) # rate of evolution trait1
)
priors_df <- data.frame(t(rbind(lapply(priors, mean), lapply(priors, sd))))
colnames(priors_df) <- c("mean", "sd")

# plot other estimated parameters
pdf(file=paste0("figures/comparison/globalparams_comparison.pdf"), height = 15, width = 18)
par(mfrow = c(11,11), mar=c(1.4,0,0,0)+0.7, mgp=c(0,0.5,0))
for(i in order(sigmasb)){
  
  est_bin <- data.frame(
    mean = as.numeric(unique(databin[[i]][["est"]][paste0(rownames(priors_df), c(".mean"))])),
    q2.5 = as.numeric(unique(databin[[i]][["est"]][paste0(rownames(priors_df), c(".q2.5"))])),
    q97.5 = as.numeric(unique(databin[[i]][["est"]][paste0(rownames(priors_df), c(".q97.5"))])))
  est_beta <- data.frame(
    mean = as.numeric(unique(databeta[[i]][["est"]][paste0(rownames(priors_df), c(".mean"))])),
    q2.5 = as.numeric(unique(databeta[[i]][["est"]][paste0(rownames(priors_df), c(".q2.5"))])),
    q97.5 = as.numeric(unique(databeta[[i]][["est"]][paste0(rownames(priors_df), c(".q97.5"))])))
  
  # scale estimates (by prior mean/std)
  for(r in 1:nrow(est_bin)){
    est_bin[r,] <- (est_bin[r,] - priors_df[r, "mean"])/priors_df[r, "sd"]
    est_beta[r,] <- (est_beta[r,] - priors_df[r, "mean"])/priors_df[r, "sd"]
  }
  rownames(est_bin) <- rownames(est_beta) <- rownames(priors_df)
  est_bin$x <- 1:6 - 0.2
  est_beta$x <- 1:6 + 0.2
  
  params <- databin[[i]][["par"]]
  
  plot.new()
  limits <- c(1-0.25,7+0.25)
  plot.window(xlim = limits, ylim = c(-3.5,3.5+.2))
  grid()
  
  labs <- c(expression(a [root]), expression(lambda  [a]), expression(sigma  [a]),
            expression(b [root]), expression(lambda  [b]), expression(sigma  [b]))
  
  param_df <- 
    data.frame(mean = c(params$a_z, params$lambda_a, params$sigma_a, params$b_z, params$lambda_b, params$sigma_b), 
               x= c(1:6))
  
  segments(y0 = param_df$mean, x0 = param_df$x-0.25, 
           y1 = param_df$mean, x1 = param_df$x+0.25,
           col = "grey50")
  segments(y0 = est_bin$q2.5, x0 = est_bin$x, 
           y1 = est_bin$q97.5, x1 = est_bin$x,
           col = "#AE2565FF")
  points(x = est_bin$x, y =  est_bin$mean, 
         pch = 19, cex = 0.5, col = "#AE2565FF")
  segments(y0 = est_beta$q2.5, x0 = est_beta$x, 
           y1 = est_beta$q97.5, x1 = est_beta$x,
           col = "#8BA749FF")
  points(x = est_beta$x, y =  est_beta$mean, 
         pch = 19, cex = 0.5, col = "#8BA749FF")
  axis(1, las = 2, cex.axis = 0.7, tck=-0.02, at = 1:6, labels = labs)
  axis(2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
  text(x = 0.7, y = 3.5+.2, label = "Binomial regression", cex = 0.45, 
       col = "#AE2565FF", adj=0)
  text(x = 0.7, y = 3+.2, label = "Beta regression", cex = 0.45, 
       col = "#8BA749FF", adj=0)
  
}
dev.off()


