
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")
palette <- colorRampPalette(kippenberger)(nspecies)

params <- list(a_z = 1, # root value intercept
               lambda_a = 0.8, # lambda intercept
               sigma_a = 0.5, # rate of evolution intercept
               b_z = 0.5, # root value trait1 slope
               lambda_b = NA, # lambda trait1
               sigma_b = NA, # rate of evolution trait1
               gamma = 5 # intercept dispersion parameter
)


pdf(file=paste0("figures/simulated_data_lambdab.sigmab.pdf"), height = 10, width = 15)
par(mfrow = c(6,11), mar=c(1,0,0,0)+1)
for(l in seq(0.0,1,0.2)){
  params$lambda_b <- l
  for(s in  seq(-0.5,0.5,0.1)){
    params$sigma_b <- s
    simulated_data <- simulate_data(spetree, params,
                                    mux = 0, sigmax = 5)
    
    plot.new()
    plot.window(xlim = c(-15,15), ylim = c(0,1))
    grid()
    points(x = simulated_data$x, y =  simulated_data$y, 
           pch = 19, cex = 0.2,
           col = palette[simulated_data$sp])
    axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    title(paste0(paste0("lambda_b=",l,", sigma_b=",s)), adj=0, cex.main = 0.7)
    
  }
  
}
dev.off()


params <- list(a_z = 1, # root value intercept
               lambda_a = NA, # lambda intercept
               sigma_a = NA, # rate of evolution intercept
               b_z = 0.5, # root value trait1 slope
               lambda_b = 0, # lambda trait1
               sigma_b = 0, # rate of evolution trait1
               gamma = 5 # intercept dispersion parameter
)

pdf(file=paste0("figures/simulated_data_lambdaa.sigmaa.pdf"), height = 10, width = 15)
par(mfrow = c(6,11), mar=c(1,0,0,0)+1)
for(l in seq(0.0,1,0.2)){
  params$lambda_a <- l
  for(s in  seq(-2,2,0.4)){
    params$sigma_a <- s
    simulated_data <- simulate_data(spetree, params,
                                    mux = 0, sigmax = 5)
    
    plot.new()
    plot.window(xlim = c(-15,15), ylim = c(0,1))
    grid()
    points(x = simulated_data$x, y =  simulated_data$y, 
           pch = 19, cex = 0.2,
           col = palette[simulated_data$sp])
    axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    title(paste0(paste0("lambda_a=",l,", sigma_a=",s)), adj=0, cex.main = 0.7)
    
  }
  
}
dev.off()

params <- list(a_z = NA, # root value intercept
               lambda_a = 0, # lambda intercept
               sigma_a = 0, # rate of evolution intercept
               b_z = NA, # root value trait1 slope
               lambda_b = 0, # lambda trait1
               sigma_b = 0, # rate of evolution trait1
               gamma = 5 # intercept dispersion parameter
)

pdf(file=paste0("figures/simulated_data_aroot.broot.pdf"), height = 15, width = 15)
par(mfrow = c(11,11), mar=c(1,0,0,0)+1)
for(a in seq(-10,10,2)){
  params$a_z <- a
  for(b in seq(-1,1,0.2)){
    params$b_z <- b
    simulated_data <- simulate_data(spetree, params,
                                    mux = 0, sigmax = 5)
    
    plot.new()
    plot.window(xlim = c(-20,20), ylim = c(0,1))
    grid()
    points(x = simulated_data$x, y =  simulated_data$y, 
           pch = 19, cex = 0.2,
           col = palette[simulated_data$sp])
    axis(1, las = 1, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    axis(2, at = c(0,0.5,1), las = 2, cex.axis = 0.7, tck=-0.02, labels=FALSE)
    title(paste0(paste0("a_z=",a,", b_z=",b)), adj=0, cex.main = 0.7)
    
  }
  
}
dev.off()




