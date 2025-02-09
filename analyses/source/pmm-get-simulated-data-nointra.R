## Code from Joly et al. 2019 ##
## He sent this in an email on 10 March 2020 ##

## This is very similar to pmm-get-simulated-data.R but the simulation removes intraspecific structure ##
## See pmm-pgls-simulations_v2.R for some of the email exchange about this code. ##

# # Parameters for troubleshooting
#nspecies=10
#nindividuals=10
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

# Simulation function (runs one rep)
one.sim.pmm.nointra <- function(nspecies, nindividuals ,B,
                        sigma.sq.x, sigma.sq.p, sigma.sq.e)
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
    spetree <- pbtree(n=nspecies,nsim=1,b=1,complete=FALSE,scale=1)
    spetree$tip.label <- letters[1:nspecies]
    # Obtain phylogenetic correlation structure (matrix)
    inter.mat <- vcv(spetree,corr=TRUE)
  
    # Make a block diagonal matrix representing the correlation structure of all species
    mat.names <- unlist(lapply(seq(1:nspecies),function(i) paste(letters[i],
                                                                 seq(1,nindividuals),sep="")))
    
    ###
    # Simulate variables
    
    # Parameters
    n.sp <- nspecies*nindividuals # number of samples
    P <- inter.mat # phylogenetic correlation matrix
    # C <- intra.mat # infraspecific correlation matrix
    # Random deviates from the normal distribution for the different effects:
    u <- rnorm(n=n.sp,mean=10) # this gets used in x (but not sure why this is set to 10 here)
    v <- rnorm(n=nspecies,mean=0) 
    # w <- rnorm(n=n.sp,mean=0) 
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
    # c <- t(chol(C*sigma.sq.c)) %*% w
    
    # residual error (not phylogenetically correlated)
    e <- matrix((sqrt(sigma.sq.e) * z), ncol=1)
    rownames(e) <- rownames(mat.names)
    
    # Calculate response variable
    y <- x %*% B + p  + e # deleted: + c
    
    
    ###
    # Data frame with the simulated data
    simdata <- list()
    simdata$data <- data.frame(y=y,x=x,animal=gsub("[0-9]*","",mat.names))
    simdata <<- simdata
    simdata$phylo <- spetree  
    simdata$intermat <- inter.mat  
    # simdata$intramat <- intra.mat  
    
  
  # Return results
  return(simdata)
}


# END
