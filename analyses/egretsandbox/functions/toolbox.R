
simulate_data_improved <- function(spetree){
  
  # parameters
  params <- list(a_z = rnorm(1,0,1.5), # root value intercept
                 lambda_a = rbeta(1,1.5,1.5), # lambda intercept
                 sigma_a = abs(rnorm(1,0,1)), # rate of evolution intercept
                 b_z = rnorm(1,0.5,1), # root value trait1 slope
                 lambda_b = rbeta(1,1.5,1.5), # lambda trait1
                 sigma_b = abs(rnorm(1,0,1)) # rate of evolution trait1
  )
  
  ## Simulate some data ## 
  # phylogenetic structure
  scaledtree_intercept <- rescale(spetree, model = "lambda", params$lambda_a)         
  params$intercepts <- fastBM(scaledtree_intercept, a = params$a_z, mu = 0, sig2 = params$sigma_a ^ 2)
  scaledtree_slope <- rescale(spetree, model = "lambda", params$lambda_b)         
  params$slopes <- fastBM(scaledtree_slope, a = params$b_z, mu = 0, sig2 = params$sigma_b ^ 2)
  
  # we consider that a study looks at only one species, but a species can be considered by several studies
  nstudies_perspecies <- round(runif(nspecies, 1, 4)) # number of different studies per species
  nstudies <- sum(nstudies_perspecies) # total number of studies
  species_study <- rep(1:nspecies, times = nstudies_perspecies) # species id considered by each study
  ntreat_perstudies <- round(runif(nstudies, 2, 6)) # number of different treatments applied per study
  nexps <- sum(ntreat_perstudies) # total number of experiments (i.e. studies*treatments)
  ntrialseeds_perexp <- round(runif(nexps, 10, 100)) # number of seeds per experiments
  studies <- rep(1:nstudies, times = ntreat_perstudies)
  species <- rep(species_study, times = ntreat_perstudies) # species id per experiments
  
  # experimental observations
  x <- rnorm(n = nexps, mean = 0, sd = 1) # treatments applied
  yhat <- plogis(params$intercepts[species] + x * params$slopes[species], location = 0, scale = 1)
  y = rbinom(n = nexps, size = ntrialseeds_perexp, prob = yhat)
  
  return(list(params = lapply(params,round,2), data = data.frame(species, studies, nseeds = ntrialseeds_perexp, x, y),
              vphy = vcv(spetree, corr = TRUE)))
  
}

get_estimates <- function(summ){
  return(
    data.frame(lambda_a.mean = summ["lambda_a" , "mean"],
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
             b_z.q97.5 = summ["b_z", "X97.5."]))
}


# based on Deirdre code
simulate_data_beta <- function(spetree, nind_species, params,
                          mux = 0, sigmax = 1,
                          location = 0, scale = 1){
  
  scaledtree_intercept <- rescale(spetree, model = "lambda", params$lambda_a)         
  intercepts <- fastBM(scaledtree_intercept, a = params$a_z, mu = 0, sig2 = params$sigma_a ^ 2)
  
  # Generate bf slope
  scaledtree_bf <- rescale(spetree, model = "lambda", params$lambda_b)         
  slopes_bf <- fastBM(scaledtree_bf, a = params$b_z, mu = 0, sig2 = params$sigma_b ^ 2)
  
  dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
  for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind_species[i]),
                       intercepts = rep(intercepts[i], nind_species[i]),
                       x = rnorm(n = nind_species[i], mux, sigmax),
                       slopes = rep(slopes_bf[i], nind_species[i])
    )
    dfhere <- rbind(dfhere, temp)
  }
  
  dfhere$yhat <- plogis(dfhere$intercepts + dfhere$x * dfhere$slopes, location = location, scale = scale)
  
  
  dfhere$shape_alpha = dfhere$yhat * params$phi
  dfhere$shape_beta = (1.0-dfhere$yhat) * params$phi
  dfhere$y = rbeta(n = nrow(dfhere), shape1 = dfhere$shape_alpha, shape2 = dfhere$shape_beta)
  dfhere$y <- pmin(pmax(dfhere$y, 1e-6), 1 - 1e-6) # small trick (yeeeek)
  
  # par(mfrow = c(1,2))
  # hist(dfhere$yhat)
  # plot(dfhere$y ~ dfhere$x, col = dfhere$sp)
  # 
  return(dfhere)

}

# based on Deirdre code
simulate_data_binomial <- function(spetree, nind_species, params, 
                                mux = 0, sigmax = 1,
                                location = 0, scale = 1){
  
  
  
  dfhere <- data.frame(sp = c(), intercept = c(), x1=numeric(), trait1=numeric())
  for (i in 1:nspecies){
    temp <- data.frame(sp = rep(i, nind_species[i]),
                       intercepts = rep(intercepts[i], nind_species[i]),
                       x = rnorm(n = nind_species[i], mux, sigmax),
                       slopes = rep(slopes_bf[i], nind_species[i])
    )
    dfhere <- rbind(dfhere, temp)
  }
  
  dfhere$yhat <- plogis(dfhere$intercepts + dfhere$x * dfhere$slopes, location = location, scale = scale)
  
  dfhere$y = rbinom(n = nrow(dfhere), size = params$ntrials, prob = dfhere$yhat)
  
  # par(mfrow = c(1,2))
  # hist(dfhere$yhat)
  # plot(dfhere$y ~ dfhere$x, col = dfhere$sp)
  # 
  return(dfhere)
  
}

plot_stuff <- function(estimates, observed, simulated_data, posteriors, comments = NULL){
  
  kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")
  palette <- colorRampPalette(kippenberger)(nspecies)
  
  palette4 <- c("#257180", "#F2E5BF", "#FD8B51", "#CB6040")
  
  plot.new()
  plot.window(xlim = c(-10,10), ylim = c(0,1))
  grid()
  points(x = simulated_data$x, y =  simulated_data$y, 
         pch = 19, cex = 0.2,
         col = palette[simulated_data$sp])
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02)
  axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.6, tck=-0.02)
  title(paste0("Observations", comments), adj=0, cex.main = 0.7)
  
  plot.new()
  limits <- c(-2,2)
  plot.window(xlim = limits, ylim = limits)
  grid()
  abline(a=0, b=1, col = "grey80")
  segments(y0 = estimates$intercepts.q2.5, x0 = observed$intercepts, 
           y1 = estimates$intercepts.q97.5, x1 = observed$intercepts,
           col = palette)
  points(x = observed$intercepts, y =  estimates$intercepts.mean, 
         pch = 19, cex = 0.5, col = palette)
  title("Intercepts", adj=0, cex.main = 0.7)
  title(ylab = "Estimated", xlab = "Observed", cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02)
  # text(simulated$intercept.mean-0.01, observed$intercepts + 0.05, 
  #      col = "grey50",
  #      labels= 1:40, cex = 0.7)
  
  plot.new()
  limits <- c(-2,2)
  plot.window(xlim = limits, ylim = limits)
  grid()
  abline(a=0, b=1, col = "grey80")
  segments(y0 = estimates$slopes.q2.5, x0 = observed$slopes, 
           y1 = estimates$slopes.q97.5, x1 = observed$slopes,
           col = palette)
  points(x = observed$slopes, y =  estimates$slopes.mean, 
         pch = 19, cex = 0.5, col = palette)
  title("Slopes", adj=0, cex.main = 0.7)
  title(ylab = "Estimated", xlab = "Observed", cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02)
  
  plot.new()
  limits <- c(posteriors$limits)
  plot.window(xlim = limits, ylim = limits)
  grid()
  points(x = posteriors[[1]], y =  posteriors[[2]], 
         pch = 19, cex = 0.2, col = palette4[posteriors$chain])
  abline(v = unique(observed$par1), lwd = 0.5)
  abline(h = unique(observed$par2), lwd = 0.5)
  title("Sampling", adj=0, cex.main = 0.7)
  title(ylab = names(posteriors)[2], xlab = names(posteriors)[1], cex.axis = 0.6, cex.lab = 0.7)
  axis(1, las = 1, cex.axis = 0.7, tck=-0.02)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02)
  
}

plot_slopes <- function(estimates, observed, params){
  
  kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")
  palette <- colorRampPalette(kippenberger)(nspecies)
  
  plot.new()
  limits <- c(-2,2)
  plot.window(xlim = limits, ylim = limits)
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
  title(paste0("a_z=", simulated_data[[i]]$params$a_z, ", l_a=", simulated_data[[i]]$params$lambda_a, ", s_a=", simulated_data[[i]]$params$sigma_a, "\n",
               "b_z=", simulated_data[[i]]$params$b_z, ", l_b=", simulated_data[[i]]$params$lambda_b, ", s_b=", simulated_data[[i]]$params$sigma_b, ", g=", simulated_data[[i]]$params$phi), 
        adj=0, cex.main = 0.5)

}





