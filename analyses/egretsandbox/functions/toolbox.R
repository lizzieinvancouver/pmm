
# based on Deirdre code
simulate_data <- function(spetree, nind_species, params,
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


