## Code from Joly et al. 2019 ##
## He sent this in an email on 10 March 2020 ##

## Lizzie asked: Do you have any simulated data + model code from the MEE paper that using just an interspecific phylogeny? Then we could do some work to better understand what the within-species centering is doing.

## He wrote: I don't have code for exactly what you want, but you can probably adapt the one I have relatively easily. In the R code attached, the function "one.sim.pmm()" simulates data and analyze it with PMM and PGLS. I think I used the same notation as in the paper. You would just have to drop the term for the intraspecific variance and simulate more values for the residual error (or add a similar term for the measurement error in addition to the residual error). Let me know if you have any other questions.

## Nacho updated this code on April 8th so it can be sourced from other pieces of code
## to generate one set of simulated data (dataset and simulated tree) according to pre-specified
## parameters (nspecies, nindividuals, strength of x-y relationship)


# # Parameters for troubleshooting
#nspecies=30
#nindividuals=40
#B = 0.25
#sigma.sq.x = 2
#sigma.sq.p = 1
#sigma.sq.c = 1
#sigma.sq.e = 1
# ms.opt="-T -I 2 5 5 -ej 0.5 2 1"


# Function to inflate a matrix (used in simulation function)
inflate.mat <- function (mat,nsp,nind) {
  inflated.mat <- matrix(apply(apply(mat,2,function(c)
    rep(c,each=nind)),1,function(c)
      rep(c,each=nind)),nrow=nsp*nind,ncol=nsp*nind)
  return(inflated.mat)
}


# Function to generate random names going beyond 27 letters
rnd.naming <- function(n = 5000) {
  a <- do.call(paste0, replicate(3, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sample(99, n, TRUE))
}


# Simulation function (runs one rep)
one.sim.pmm <- function(nspecies, nindividuals,B, scaletree,
                        sigma.sq.x,sigma.sq.p,sigma.sq.c,sigma.sq.e,
                        phyloSlope)
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
  
   
  # Simulate species tree with a pure birth model
    spetree <- pbtree(n=nspecies,nsim=1,b=1,complete=FALSE,scale=scaletree)
    spetree$tip.label <- rnd.naming(nspecies)
    # Obtain phylogenetic correlation structure (matrix)
    inter.mat <- vcv(spetree, corr=TRUE)
  
   
    # we generate here a list of star-like genealogies (vcvs)
    # since we are not concerned (so far) with the genetic relationship among individuals
    matempty <- matrix(rep(0,nindividuals*nindividuals),ncol=nindividuals,nrow=nindividuals,)
    diag(matempty) <- rep(1, nindividuals)
    
    popstruct <- lapply(seq(1:nspecies), function(c)c=matempty)
    # popstruct <- lapply(seq(1:nspecies),function(c) matrix(apply(sapply(sim.ms(nsam=nindividuals,nreps=ngen,ms.command=ms.opt),function(c) {xx <- vcv(c,corr=TRUE);or<-order(as.numeric(gsub("[a-z]","",colnames(xx))));return(xx[or,or])}),1,mean),ncol=nindividuals,nrow=nindividuals))

    # Make a block diagonal matrix representing the correlation structure of all species

    mat.names <- unlist(lapply(seq(1:nspecies),function(i) paste(spetree$tip.label[i],
                                                                 seq(1,nindividuals),sep=".")))
    
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

    # What if we want the phylogeny to structure the slopes? (By Lizzie)
    if(phyloSlope){
    B1 <- B 
    B <-  t(chol(P*B)) %*% v
    B <- matrix(rep(B[,1],times=rep(nindividuals,nspecies)), ncol=1)
    rownames(B) <- mat.names

    # Calculate response variable
    # Sidenote that identical(x%*%B1, x*B1)
    y <- x * B + p + c + e
    }
    if(!phyloSlope){
    # Calculate response variable
    y <- x %*% B + p + c + e
    }
    
    
    ###
    # Data frame with the simulated data
    simdata <- list()
    simdata$data <- data.frame(y=y,x=x,animal=unlist(lapply(strsplit(mat.names,"\\.",fixed=F),function(x)x[1])))
    simdata <<- simdata
    simdata$phylo <- spetree  
    simdata$intermat <- inter.mat  
    simdata$intramat <- intra.mat  
    
  
  # Return results
  return(simdata)
}


# END
