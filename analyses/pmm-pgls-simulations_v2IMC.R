## Code from Joly et al. 2019 ##
## He sent this in an email on 10 March 2020 ##

## Lizzie asked: Do you have any simulated data + model code from the MEE paper that using just an interspecific phylogeny? Then we could do some work to better understand what the within-species centering is doing.

## He wrote: I don't have code for exactly what you want, but you can probably adapt the one I have relatively easily. In the R code attached, the function "one.sim.pmm()" simulates data and analyze it with PMM and PGLS. I think I used the same notation as in the paper. You would just have to drop the term for the intraspecific variance and simulate more values for the residual error (or add a similar term for the measurement error in addition to the residual error). Let me know if you have any other questions.



# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} 

## load packages
library(shinystan)
library(caper)
library(brms)
library(pez)
library(rstan)
library(phytools)
library(MCMCglmm)
library(dplyr)
library(knitr)
library(broom)
library(loo)
# Packages required
require(ape)
require(nlme)
require(Matrix)
require(phyclust)

source("corIntra.R")

# # Parameters for troubleshooting
nspecies=10
nindividuals=10
ngen=20
B = 0.25
sigma.sq.x = 2
sigma.sq.p = 1
sigma.sq.c = 1
sigma.sq.e = 1
# ms.opt="-T -I 2 5 5 -ej 0.5 2 1"


# Function to inflate a matrix (used in simulation function)
inflate.mat <- function (mat,nsp,nind) {
  inflated.mat <- matrix(apply(apply(mat,2,function(c)
    rep(c,each=nind)),1,function(c)
      rep(c,each=nind)),nrow=nsp*nind,ncol=nsp*nind)
  return(inflated.mat)
}

# Simulation function (runs one rep)
one.sim.pmm <- function(nspecies=10,nindividuals=10,ngen=20,B = 0.25,
                        sigma.sq.x = 2,sigma.sq.p = 1,sigma.sq.c = 1,sigma.sq.e = 1,
                        ms.opt="-T -I 2 5 5 -ej 0.5 2 1")
{
  
  ## Global parameters
  #
  # nspecies: number of species in the tree
  # nindividuals: number of individuals sampled per species
  # ngen: number of gene genealogies from which to estimate the infraspecific
  #         genetic structure.
  # ms.opt: Command line options for running the 'ms' coalescent simulation 
  #           program
  #
  ## Character simulation parameters
  #
  # sigma.sq.x: Rate of evolution in x
  # sigma.sq.p: Phylogenetic variance
  # sigma.sq.c: Intraspecific correlation structure variance
  # sigma.sq.e: Rate of residual variation
  # B: Regression slope
  #
  
  # Packages required
  require(ape)
  require(nlme)
  require(Matrix)
  require(phytools)
  require(MCMCglmm)
  require(phyclust)
  
  # To check for PGLS convergence
  convergence=FALSE
  times=0
  while (convergence==FALSE) {
    # Simulate species tree with a pure birth model
    spetree <- pbtree(n=nspecies,nsim=1,b=1,complete=FALSE,scale=1)
    spetree$tip.label <- letters[1:nspecies]
    # Obtain phylogenetic correlation structure (matrix)
    inter.mat <- vcv(spetree,corr=TRUE)
  
    # Function to perform coalescent simulations using the 'ms' software
    sim.ms <- function(nsam=nindividuals,nreps=20,ms.command="-T") {
      require(phyclust)
      ret.ms <- ms(nsam = nsam, nreps = nreps, opts = ms.command)
      return(read.tree(text = ret.ms,skip=1,comment="/"))
    }
  
    # Simulate an infraspecific correlation structure for each species. The 
    # correlation structure is the mean of the VCV matrices of 20 gene genealogies
    #popstruct <- lapply(seq(1:nspecies),function(c) matrix(apply(sapply(sim.ms(nsam=nindividuals,nreps=ngen,ms.command=ms.opt),
    #                                                                    function(c) {xx <- vcv(c,corr=TRUE);
     #                                                                   or<-order(as.numeric(gsub("[a-z]","",colnames(xx))));
      #                                                                  return(xx[or,or])})
       #                                                          ,1,mean),
        #                                                   ncol=nindividuals,nrow=nindividuals))
    
    # we can substutite the previous code to get a list of star-like genealogies (vcvs)
    # since we do not care about the genetic relationship among individuals
    matempty <- matrix(rep(0,nspecies^2),ncol=nspecies,nrow=nspecies,)
    diag(matempty) <- rep(1, nspecies)
    
    popstruct <- lapply(seq(1:nspecies), function(c)c=matempty)
     
    
    
    
    # Make a block diagonal matrix representing the correlation structure of all species
    mat.names <- unlist(lapply(seq(1:nspecies),function(i) paste(letters[i],
                                                                 seq(1,nindividuals),sep="")))
    intra.mat <- matrix(bdiag(lapply(seq(1,length(popstruct)),function(c) 
      popstruct[[c]])),nrow=(nspecies*nindividuals),dimnames=list(mat.names,mat.names))
    
    
    ###
    # Simulate variables
    
    # Parameters
    n.sp <- nspecies*nindividuals # number of samples
    P <- inter.mat # phylogenetic correlation matrix
    C <- intra.mat # infraspecific correlation matrix
    # Random deviates from the normal distribution for the different effects:
    u <- rnorm(n=n.sp,mean=10) 
    v <- rnorm(n=nspecies,mean=0) 
    w <- rnorm(n=n.sp,mean=0) 
    z <- rnorm(n=n.sp,mean=0) 
    
    # Get values for x (not phylogenetically correlated)
    x <- matrix((sqrt(sigma.sq.x) * u), ncol=1)
    rownames(x) <- mat.names
    
    # Phylogenetic random effect
    p <- t(chol(P*sigma.sq.p)) %*% v
    # Expand for all individuals within species
    p <- matrix(rep(p[,1],times=rep(nindividuals,nspecies)),ncol=1)
    rownames(p) <- mat.names
    
    # Intraspecific random effect
    c <- t(chol(C*sigma.sq.c)) %*% w
    
    # residual error (not phylogenetically correlated)
    e <- matrix((sqrt(sigma.sq.e) * z), ncol=1)
    rownames(e) <- rownames(mat.names)
    
    # Calculate response variable
    y <- x %*% B + p + c + e
    
    
    ###
    # Data frame with the simulated data
    simdata <- list()
    simdata$data <- data.frame(y=y,x=x,animal=gsub("[0-9]*","",mat.names))
    simdata <<- simdata
  
    ###
    # Prepare infraspecific correlation matrix for MCMCglmm
    #
    # Single value decomposition of the intraspecific structure
    #   -> follow code of Stone et al. 2012
    intra.svd <- svd(intra.mat)
    intra.svd <- intra.svd$v %*% (t(intra.svd$u) * sqrt(intra.svd$d))
    rownames(intra.svd) <- colnames(intra.svd) <- rownames(intra.mat)
    # The following is important and MCMCglmm searches the variables in the 
    # global environment
    intra.svd <<- intra.svd
    
    ###
    # Run the MCMCpglmm analyses with different correlation structures.
    # For each, we will use a very diffuse prior.
    
    # Model M.0 -> no random effects
    priorpr.m0 <- list(R = list(V = 1, nu = 0.002))
    M.0 <- MCMCglmm(y ~ x, 
                    data=simdata$data, scale=TRUE,prior=priorpr.m0,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
    
    # Model M.1 -> only phylogenetic random effects
    priorpr.m1 <- list(R = list(V = 1, nu = 0.002), 
                       G = list(G1 = list(V = 1, nu = 0.002)))
    M.1 <- MCMCglmm(y ~ x, random = ~ animal, pedigree = spetree,
                    data=simdata$data, scale=TRUE,prior=priorpr.m1,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
  
    # Model M.2 -> only infraspecific random effects
    M.2 <- MCMCglmm(y ~ x, random =~ idv(intra.svd), 
                    data=simdata$data, scale=TRUE,prior=priorpr.m1,
                    nitt=21000,thin=10,burnin=1000,verbose=FALSE)
    
    # Model M.3 -> infraspecific and interspecific random effects
    priorpr.m3 <- list(R = list(V = 1, nu = 0.002), 
                        G = list(G1 = list(V = 1, nu = 0.002),
                                 G2 = list(V = 1, nu = 0.002)))
    M.3 <- MCMCglmm(y ~ x, random =~ animal + idv(intra.svd), pedigree = spetree,
                     data=simdata$data, scale=TRUE,prior=priorpr.m3,
                     nitt=21000,thin=10,burnin=1000,verbose=FALSE)
  
    ###
    # Run PGLS
    
    inter.mat.inflated <- inflate.mat(inter.mat,nspecies,nindividuals)
    # Phylogenetic correlation matrices
    #
    # The intra correlation structure takes the inter-specific and the intra-
    # specific correlation matrices as input.
    pgls.intra.cor <- corIntra(0, inter=inter.mat.inflated, intra=intra.mat, fixed=TRUE)
    pgls.inter.cor <- corIntra(1, inter=inter.mat.inflated, intra=intra.mat, fixed=TRUE)
    pgls.intrainter.cor <- corIntra(0.5, inter=inter.mat.inflated, intra=intra.mat, fixed=FALSE)
    
    # Fit PGLS
    gls.mod0 <- gls.mod1 <- gls.mod2 <- gls.mod3 <-NULL
    try(gls.mod0 <- gls(y~x,data=simdata$data)) # no phylogenetic correlation
    try(gls.mod1 <- gls(y~x,correlation = pgls.inter.cor, data=simdata$data),silent=TRUE) # inter structure 
    try(gls.mod2 <- gls(y~x,correlation = pgls.intra.cor, data=simdata$data),silent=TRUE) # intra structure
    try(gls.mod3 <- gls(y~x,correlation = pgls.intrainter.cor, data=simdata$data),silent=TRUE) # Estimating delta 
    
    if(any(is.null(gls.mod0),is.null(gls.mod1),is.null(gls.mod2),is.null(gls.mod3))) {
      times=times+1
      if (times>=10) return(NULL)
    } else {
      convergence=TRUE
    }
  }
  
  ###
  # Compile the results
  
  # Model summaries
  M.0.sum <- summary(M.0)
  M.1.sum <- summary(M.1)
  M.2.sum <- summary(M.2)
  M.3.sum <- summary(M.3)
  gls.mod0.sum <- summary(gls.mod0)
  gls.mod1.sum <- summary(gls.mod1)
  gls.mod2.sum <- summary(gls.mod2)
  gls.mod3.sum <- summary(gls.mod3)
  
  # Data frame
  res <- data.frame(models=c("M.0","M.1","M.2","M.3","gls0","gls1","gls2","gls3"),
                    #random.effects=c("NA","inter","infra","inter+infra","ols","glsInter","glsIntra","glsInterIntra"),
                    DIC=c(M.0$DIC,M.1$DIC,M.2$DIC,M.3$DIC,-999,-999,-999,-999),
                    slope=c(M.0.sum$solutions[2,1],M.1.sum$solutions[2,1],
                            M.2.sum$solutions[2,1],M.3.sum$solutions[2,1],
                            gls.mod0$coefficients[2],gls.mod1$coefficients[2],
                            gls.mod2$coefficients[2],gls.mod3$coefficients[2]),
                    slope.sd=c(sd(M.0$Sol[,'x']),sd(M.1$Sol[,'x']),
                               sd(M.2$Sol[,'x']),sd(M.3$Sol[,'x']),
                               gls.mod0.sum$tTable[2,2],
                               gls.mod1.sum$tTable[2,2],
                               gls.mod2.sum$tTable[2,2],
                               gls.mod3.sum$tTable[2,2]),
                    slope.signif=c(ifelse(M.0.sum$solution[2,5]<0.05,1,0),
                                   ifelse(M.1.sum$solution[2,5]<0.05,1,0),
                                   ifelse(M.2.sum$solution[2,5]<0.05,1,0),
                                   ifelse(M.3.sum$solution[2,5]<0.05,1,0),
                                   ifelse(gls.mod0.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod1.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod2.sum$tTable[2,4]<0.05,1,0),
                                   ifelse(gls.mod3.sum$tTable[2,4]<0.05,1,0)))
  
  # Calculate the absolute difference with the 'true' slope
  res$slope.diff <- abs(B-res$slope)
  res$best.slope <- ifelse(res$slope.diff==min(res$slope.diff),1,0)
  
  # Calculate the variance contributions
  var.contr <- function(model){
    xx <- model$VCV/apply(model$VCV,1,sum)
    apply(xx,2,mean)
  }
  M.0.var <- var.contr(M.0)
  M.1.var <- var.contr(M.1)
  M.2.var <- var.contr(M.2)
  M.3.var <- var.contr(M.3)
  
  # Proportion of intraspecific structure relative to the complete genetic structure
  res$intra.prop <- c(NA,NA,NA,M.3.var["intra.svd."]/(M.3.var["intra.svd."]+M.3.var["animal"]),
                      NA,NA,NA,(1-gls.mod3.sum$modelStruct$corStruct[1]))
  
  # Heritability (h^2)
  res$h2 <- c(NA,M.1.var["animal"]/sum(M.1.var),
              M.2.var["intra.svd."]/sum(M.2.var),
              (M.3.var["intra.svd."]+M.3.var["animal"])/sum(M.3.var),
              NA,NA,NA,NA)
  
  ###
  # Return results
  return(res)
}


###
# Simulations for a whole set of settings

library(parallel)
# Parrameters to test
B.vals <- c(0,0.1,0.25)
sigma.sq.b.vals <- c(0.5, 1, 1.5)
#B.vals <- c(0.1)
#sigma.sq.b.vals <- c(0.5, 1.5)
nsim=500
no_cores <- detectCores() - 1 # Calculate the number of cores

#test (sapply)
#results <- sapply(c(1:2),function(x) one.sim.pmm(nspecies=10,nindividuals=10,B=0.1,ngen=20,sigma.sq.x=2,sigma.sq.p=1,sigma.sq.c=1,sigma.sq.e=1,ms.opt="-T -I 3 3 3 4 -ej 0.3 3 2 -ej 0.8 2 1"))

# table to place the results
res.table <- data.frame(model=character(),
                        beta=numeric(),sigma.a=numeric(),sigma.b=numeric(),
                        slope=numeric(),slope.sd=numeric(),accuracy=numeric(),
                        precision=numeric(),power=numeric(),
                        intra.prop=numeric(),intra.prop.sd=numeric(),
                        h2=numeric(),h2.sd=numeric(),
                        best=numeric())

for (B.val in B.vals) {
  for (sigma.sq.b.val in sigma.sq.b.vals) {
    # Run simulations
    cl <- makeCluster(no_cores,type="FORK",outfile="sim.txt") # Initiate cluster
    results <- parSapply(cl,c(1:nsim),function(x) one.sim.pmm(nspecies=10,nindividuals=10,B=B.val,ngen=20,sigma.sq.x=2,sigma.sq.p=(2-sigma.sq.b.val),sigma.sq.c=sigma.sq.b.val,sigma.sq.e=1,ms.opt="-T -I 3 3 3 4 -ej 0.3 3 2 -ej 0.8 2 1"))
    stopCluster(cl)
    # Extract stats
    slope.mean <- apply(array(unlist(results["slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    slope.sd <- apply(array(unlist(results["slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    slope.diff <- apply(array(unlist(results["slope.diff",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    slope.prec <- apply(array(unlist(results["slope.sd",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    prop.best <- apply(array(unlist(results["best.slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    prop.signif <- apply(array(unlist(results["slope.signif",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    intra.mean <- apply(array(unlist(results["intra.prop",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    intra.sd <- apply(array(unlist(results["intra.prop",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    h2.mean <- apply(array(unlist(results["h2",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    h2.sd <- apply(array(unlist(results["h2",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    dic <- apply(array(unlist(results["DIC",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    # Store values in temp dataframe
    res.temp <- data.frame(model=c("M0","M1","M2","M3","gls0","gls1","gls2","gls3"),
                           beta=B.val,sigma.a=(2-sigma.sq.b.val),sigma.b=sigma.sq.b.val,
                           slope=slope.mean,slope.sd=slope.sd,accuracy=slope.diff,
                           precision=slope.prec, power=prop.signif,
                           intra.prop=intra.mean,intra.prop.sd=intra.sd,
                           h2=h2.mean,h2.sd=h2.sd,
                           best=c(prop.best[1],prop.best[2],prop.best[3],prop.best[4],prop.best[5],prop.best[6],prop.best[7],prop.best[8]))
    # Paste values in the table of results
    res.table<-rbind(res.table,res.temp)
    cat("  -> Done sigma.b =",sigma.sq.b.val,"\n")
  }
  cat("Done Beta =",B.val,"\n")
}

### With 20 species and 20 ind per species
for (B.val in B.vals) {
  for (sigma.sq.b.val in sigma.sq.b.vals) {
    # Run simulations
    cl <- makeCluster(no_cores,type="FORK",outfile="sim.txt") # Initiate cluster
    results <- parSapply(cl,c(1:nsim),function(x) one.sim.pmm(nspecies=20,nindividuals=20,B=B.val,ngen=20,sigma.sq.x=2,sigma.sq.p=(2-sigma.sq.b.val),sigma.sq.c=sigma.sq.b.val,sigma.sq.e=1,ms.opt="-T -I 3 6 7 7 -ej 0.3 3 2 -ej 0.8 2 1"))
    stopCluster(cl)
    # Extract stats
    slope.mean <- apply(array(unlist(results["slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    slope.sd <- apply(array(unlist(results["slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    slope.diff <- apply(array(unlist(results["slope.diff",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    slope.prec <- apply(array(unlist(results["slope.sd",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    prop.best <- apply(array(unlist(results["best.slope",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    prop.signif <- apply(array(unlist(results["slope.signif",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    intra.mean <- apply(array(unlist(results["intra.prop",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    intra.sd <- apply(array(unlist(results["intra.prop",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    h2.mean <- apply(array(unlist(results["h2",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    h2.sd <- apply(array(unlist(results["h2",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,sd)
    dic <- apply(array(unlist(results["DIC",]),dim=c(8,dim(results)[2]),dimnames=list(as.vector(results[1,1][[1]]),NULL)),1,mean)
    # Store values in temp dataframe
    res.temp <- data.frame(model=c("M0","M1","M2","M3","gls0","gls1","gls2","gls3"),
                           beta=B.val,sigma.a=(2-sigma.sq.b.val),sigma.b=sigma.sq.b.val,
                           slope=slope.mean,slope.sd=slope.sd,accuracy=slope.diff,
                           precision=slope.prec, power=prop.signif,
                           intra.prop=intra.mean,intra.prop.sd=intra.sd,
                           h2=h2.mean,h2.sd=h2.sd,
                           best=c(prop.best[1],prop.best[2],prop.best[3],prop.best[4],prop.best[5],prop.best[6],prop.best[7],prop.best[8]))
    # Paste values in the table of results
    res.table<-rbind(res.table,res.temp)
    cat("  -> Done sigma.b =",sigma.sq.b.val,"\n")
  }
  cat("Done Beta =",B.val,"\n")
}

save(res.table,file="results_pgls2.Rdata")


# Summarize simulations results

res.table
load(file="results_pgls2.Rdata")

levels(res.table$model) <- c("OLS","PGLS inter","PGLS intra","PGLS inter + intra","PMM null","PMM inter","PMM intra","PMM inter + intra")

###
# Figures

library(ggplot2)

# power
fig.pow <- ggplot(res.table,aes(y=power,x=factor(paste(sigma.a,sigma.b,sep=" : ")), group=model)) +
  geom_line(aes(colour=model)) +
  geom_point(aes(shape=model, colour=model, fill=model),size=2) +
  xlab(expression(paste(sigma[a]," : ",sigma[b]))) +
  ylab("Proportion of significant results") +
  facet_wrap(~ beta, labeller = label_bquote(beta == .(beta))) + theme_light() +
  theme(strip.background=element_rect(fill="gray85"), strip.text=element_text(size=10,colour="black")) +
  scale_shape_manual(values = c(2,5,0,1,24,23,22,21)) +
  scale_fill_brewer(type = "qualitative", palette = "Dark2", direction = -1) +
  scale_colour_brewer(type = "qualitative", palette = "Dark2", direction = -1)
fig.pow

# accuracy
fig.acc <- ggplot(subset(res.table,beta == 0.25),aes(y=accuracy,x=factor(paste(sigma.a,":",sigma.b)), group=model)) +
  geom_line(aes(colour=model)) +
  geom_point(aes(shape=model, colour=model, fill=model),size=2) +
  ylab("") +
  xlab(expression(paste(sigma[a]," : ",sigma[b]))) +
  facet_wrap(~ beta, labeller = label_bquote(Accuracy)) + theme_light() + 
  theme(legend.position="none", strip.background=element_rect(fill="gray85"), strip.text=element_text(size=10,colour="black")) +
  scale_colour_brewer(type = "qualitative", palette = "Dark2", direction = -1) +
  scale_shape_manual(values = c(2,5,0,1,24,23,22,21)) +
  scale_fill_brewer(type = "qualitative", palette = "Dark2", direction = -1)
fig.acc

# precision
fig.pre <- ggplot(subset(res.table,beta == 0.25),aes(y=precision,x=factor(paste(sigma.a,sigma.b,sep=" : ")), group=model)) +
  geom_line(aes(colour=model)) +
  geom_point(aes(shape=model, colour=model, fill=model),size=2) +
  ylab("") +
  xlab(expression(paste(sigma[a]," : ",sigma[b]))) +
  facet_wrap(~ beta, labeller = label_bquote(Precision)) + theme_light() + 
  theme(strip.background=element_rect(fill="gray85"), strip.text=element_text(size=10,colour="black")) +
  scale_shape_manual(values = c(2,5,0,1,24,23,22,21)) +
  scale_fill_brewer(type = "qualitative", palette = "Dark2", direction = -1) +
  scale_colour_brewer(type = "qualitative", palette = "Dark2", direction = -1)
fig.pre

# heredity
fig.h2 <- ggplot(subset(res.table,beta == 0.25),aes(y=h2,x=factor(paste(sigma.a,sigma.b,sep=" : ")), group=model)) +
  geom_hline(aes(yintercept = 0.666),lty=3) +
  geom_line(aes(colour=model)) +
  geom_point(aes(shape=model, colour=model, fill=model),size=2) +
  #ggtitle("C) Heredity") + 
  #ylab(expression(paste("heredity (",h^2,")"))) +
  ylab("") +
  xlab(expression(paste(sigma[a]," : ",sigma[b]))) +
  facet_wrap(~ beta, labeller = label_bquote(Heredity)) + theme_light() +
  theme(strip.background=element_rect(fill="gray85"), strip.text=element_text(size=10,colour="black")) +
  scale_shape_manual(values = c(2,5,0,1,24,23,22,21)) +
  scale_fill_brewer(type = "qualitative", palette = "Dark2", direction = -1) +
  scale_colour_brewer(type = "qualitative", palette = "Dark2", direction = -1)
fig.h2

library(grid)
pdf("simuls_pgls.pdf",width = 6, height = 3)
grid.draw(cbind(ggplotGrob(fig.acc), ggplotGrob(fig.pre), size="last"))
dev.off()

pdf("simuls_power_pgls.pdf",width = 8, height = 3)
fig.pow
dev.off()

# END
