

smbin_mod <- stan_model("/home/victor/projects/phylo/binomiallikelihood_expquad.stan")

simulated <- list(
  params = lapply(params,round,2), 
  data = data.frame(species, studies, nseeds = ntrialseeds_perexp, y),
  vphy = round(cophenetic(spetree), 3))
  # vphy = round(maxlength-vcv(spetree, corr = FALSE),3)) # 

mdl.data <- list(y = simulated$data$y,
                 N = nrow(simulated$data),
                 Nsp = nspecies,
                 sp = simulated$data$species,
                 Ntrials = simulated$data$nseeds,
                 dist= simulated$vphy)

fit_expquad <- sampling(smbin_mod, mdl.data, 
                iter = 4000, warmup = 3000,
                chains = 4, cores = 4)


summ_expquad <- data.frame(summary(fit_expquad)[["summary"]])
summ <- data.frame(summary(fit)[["summary"]])

exp_quad <- function(dx, gamma, rho){
  return(gamma^2 * exp(-1/2*(dx/rho)^2))
}



plotdf <- data.frame(dist = seq(0, max(cophenetic(spetree)), 0.1), cov = NA)
plotdf$cov_expquad <- maxlength*sapply(plotdf$dist, function(i) exp_quad(i, gamma = summ_expquad["sigma_a","mean"], rho = summ_expquad["rho","mean"]))



par(mfrow = c(1,1))
plot.new()
plot.window(xlim = c(0,max(distance)), ylim = c(0,maxlength))
grid()
axis(1, las = 1, cex.axis = 1, tck=-0.02, at = c(0,max(distance)), labels=c("0", expression(paste("2", varphi))))
axis(2, las = 1, cex.axis = 1, tck=-0.02, at = c(0,params$lambda_a*maxlength,maxlength), labels=c("0", expression(paste(lambda, "*", varphi)), expression(paste(varphi))))
title(main = "Correlogram", xlab = expression(paste(Delta,x)),  line=1, cex.main = 1)
lines(plotdf$dist, plotdf$cov_expquad, col="white", lwd=4)
lines(plotdf$dist, plotdf$cov_expquad, col="#6E8537FF", lwd=2, lty = 2)
lines(dataplot_sub$dist, dataplot_sub$cov, col="white", lwd=4)
lines(dataplot_sub$dist, dataplot_sub$cov, col="#AE2565FF", lwd=2)
points(dataplot[dataplot$dist == 0, "dist"], dataplot[dataplot$dist == 0, "cov"], 
       col="#AE2565FF", pch = 20)






intercepts.mean = summ[paste0("a[", 1:nspecies, "]"), "mean"]
intercepts.mean_expquad = summ_expquad[paste0("a[", 1:nspecies, "]"), "mean"]

par(mfrow = c(1,2))
plot(intercepts.mean_expquad ~ params$intercepts, xlab = "True data", ylab = "Exp. quad. model")
abline(0,1)
plot(intercepts.mean ~ params$intercepts, xlab = "True data", ylab = "True generating process")
abline(0,1)

