## Started 18 April 2020 ##
## By Lizzie to start ##
## Built off testdata_generate.R in OSPREE repo ##


set.seed(73)

library(msm) # for truncated normal distribution

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Ignacio", getwd())>0)) { 
  setwd("~/GitHub/pmm/analyses/") 
} else
setwd("~/Documents/git/teaching/stan/pmm/analyses")


##################################################################
## This version has only simple linear model:
# bb ~ force
# species varies by slopes and intercepts (aka random slopes and intercepts)
# centered version (all predictors are centered at 0)
##################################################################

# Code designed with normal distribution for intercepts for a and sigma_y

# nlab = 10 # number of labgroups
nsp = 30 # number of species
ntot = 50 # numbers of obs per species. 

#  with species  (note to self: This is not the best, better to draw from a distribution)

#  with species  (note to self: This is not the best, better to draw from a distribution)
intermean <- 30 # mean for selecting intercept (days to BB) across all species
intersd <- 3 # SD for selecting species intercepts
spint <- rnorm(nsp, intermean, intersd)  # different intercepts by species

# now start building ...
testdat2 <- vector()

# assumptions:
# (a) predictors are centered
# (b) predictors are not correlated
# (d) each obs is a different set of treatments

# and some important points ...
# (z) the below draws treatments from distribution in such a way that there is a lot more variation than we have


for(i in 1:nsp){ # loop over species. i = 1

    # continuous predictors, generate level (if you will) for each observation
    force = rnorm(ntot, 0, 2)
    
    # set up effect sizes
    forcecoef = -2

    # SD for each treatment
    forcecoef.sd = 0.5 

    # build model matrix 
    mm <- model.matrix(~force, data.frame(force))

    # coefficients need to match the order of the colums in the model matrix (mm)
    # so here, that's intercept, chill, force, photo
    coeff <- c(spint[i], 
             rnorm(1, forcecoef, forcecoef.sd)
             )
  
    bb <- rnorm(n = ntot, mean = mm %*% coeff, sd = 0.1)
  
    testdatx2 <- data.frame(bb, sp = i, force)
  
    testdat2 <- rbind(testdat2, testdatx2)  
}

#summary(lm(bb ~ force, data = testdat2)) # sanity check

#library(lme4)
#summary(lmer(bb ~ force + (1|sp), data = testdat2))


##
# try the model 
datalist.bb <- with(testdat2, 
    list(y = bb, # bbint 
         force = as.numeric(force), 
         sp = as.numeric(sp),
         N = nrow(testdat2),
         n_sp = length(unique(sp))
         )
)
 
osp.bb <- stan('stan/nointer_2level_force.stan', data = datalist.bb, 
                 iter = 2000, chains=4
                  ) # seems good!
