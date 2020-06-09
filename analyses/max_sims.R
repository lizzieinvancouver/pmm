# PGLMM Simulations

require(phytools)
require(geiger)

# Exploring vcv calculations
tree <- pbtree(n=10)

# Approach from Statistical Rethinking 2nd Ed ()
Rbm <- corBrownian(phy=tree)
V <- vcv(Rbm)
V

# convert to correlation matrix
R <- V / max(V)
R


# VCV for star phylogenies
star_tree <- rescale(tree, model="lambda", lambda=0)
plot(star_tree)
vcv(star_tree, cor=T)


# Comparing vcv() to vcv(, cor=TRUE)
require(brms)

# Simulating data
n_sp <- 100
tree <- pbtree(n=n_sp)
plot(tree)
V <- vcv(tree, cor=T)


sim.trait1 <- sim.char(tree, 1, root=0, model="BM")[,,1]
non_phylo_error1 <- rnorm(n_sp, mean=0, sd=0)

sim.trait1 <- sim.trait1 + non_phylo_error1
sd_fit1 <- MASS::fitdistr(sim.trait1, "normal")$estimate[2]


# sim.trait2 <- sim.char(tree, 1, root=0, model="BM")[,,1]
sim.trait2 <- sim.trait1 * 2.2
MASS::fitdistr(sim.trait2, "normal")

non_phylo_error2 <- rnorm(n_sp, mean=0, sd=sd_fit1)
sim.trait2 <- sim.trait2 + non_phylo_error2
MASS::fitdistr(sim.trait2, "normal")

# plot(sim.trait ~ sim.trait2)

dat <- data.frame(species=tree$tip.label, trait1=sim.trait1, trait2=sim.trait2)


# PGLS
require(caper)

dat.caper <- comparative.data(phy=tree, data=dat, names.col="species")
pgls1 <- pgls(trait2 ~ trait1, data=dat.caper, lambda="ML")
summary(pgls1)


# BRMS

m1 <- brm(trait2 ~ trait1 + (1|species),
    prior = c(prior(normal(0, 3), class = Intercept),
                prior(normal(0, 3), class = b),
                prior(normal(0, 3), class = sd)),
  data=dat, family="gaussian", cov_ranef = list(species = V),
  iter=3000, cores=4,
  control=list(adapt_delta=0.95, max_treedepth=12))

m1
pairs(m1)

# phylogenetic heritability 
hyp <- paste("sd_species__Intercept^2 /", "(sd_species__Intercept^2 + sigma) = 0")
hyp1 <- hypothesis(m1, hyp, class = NULL)

plot(density(hyp1$samples$H1))

lambda_1 <- pgls.profile(pgls1)
plot(lambda_1)




# Trying to simulate data with a known lambda value

# Simulating data
n_sp <- 100
tree <- pbtree(n=n_sp)
plot(tree)
V <- vcv(tree, cor=T)

lambda_tree <- rescale(tree, model="lambda", 0.5)

sim.trait1 <- sim.char(lambda_tree, 1, root=0, model="BM")[,,1]
non_phylo_error1 <- rnorm(n_sp, mean=0, sd=0)

sim.trait1 <- sim.trait1 + non_phylo_error1
sd_fit1 <- MASS::fitdistr(sim.trait1, "normal")$estimate[2]


sim.trait2 <- sim.char(lambda_tree, 1, root=0, model="BM")[,,1]
# sim.trait2 <- sim.trait1 * 1
# MASS::fitdistr(sim.trait2, "normal")

non_phylo_error2 <- rnorm(n_sp, mean=0, sd=0)
sim.trait2 <- sim.trait2 + non_phylo_error2
MASS::fitdistr(sim.trait2, "normal")

# plot(sim.trait2 ~ sim.trait1)
# cor(sim.trait2, sim.trait1)

dat <- data.frame(species=tree$tip.label, trait1=sim.trait1, trait2=sim.trait2)


# PGLS
require(caper)

# Fitting with the original tree
dat.caper <- comparative.data(phy=tree, data=dat, names.col="species")
pgls1 <- pgls(trait2 ~ trait1, data=dat.caper, lambda="ML")
summary(pgls1)
plot(pgls.profile(pgls1))

# Fitting with the lambda tree
dat.caper.lt <- comparative.data(phy=lambda_tree, data=dat, names.col="species")
pgls.lt <- pgls(trait2 ~ trait1, data=dat.caper.lt, lambda="ML")
summary(pgls.lt)
plot(pgls.profile(pgls.lt))


# Looks like we can recapture lambda when two intependent traits are simulated from it

# BRMS
m2 <- brm(trait2 ~ trait1 + (1|species),
    prior = c(prior(normal(0, 3), class = Intercept),
                prior(normal(0, 3), class = b),
                prior(normal(0, 3), class = sd)),
  data=dat, family="gaussian", cov_ranef = list(species = V),
  iter=3000, cores=4,
  control=list(adapt_delta=0.95, max_treedepth=12))

m2
pairs(m2)

# phylogenetic heritability 
hyp <- paste("sd_species__Intercept^2 /", "(sd_species__Intercept^2 + sigma) = 0")
hyp2 <- hypothesis(m2, hyp, class = NULL)

str(hyp2)

plot(density(hyp2$samples$H1), xlim=c(0,1))
abline(v=hyp2$hypothesis$Estimate, col=4, lty=2)
abline(v=hyp2$hypothesis$CI.Lower, col=2, lty=2)
abline(v=hyp2$hypothesis$CI.Upper, col=2, lty=2)

lambda_1 <- pgls.profile(pgls1)
plot(lambda_1)


# Trying with VCV cor = FALSE
V <- vcv(tree, cor=F)

# BRMS
m3 <- brm(trait2 ~ trait1 + (1|species),
    prior = c(prior(normal(0, 3), class = Intercept),
                prior(normal(0, 3), class = b),
                prior(normal(0, 3), class = sd)),
  data=dat, family="gaussian", cov_ranef = list(species = V),
  iter=3000, cores=4,
  control=list(adapt_delta=0.95, max_treedepth=12))

m3
pairs(m3)

# phylogenetic heritability 
hyp <- paste("sd_species__Intercept^2 /", "(sd_species__Intercept^2 + sigma) = 0")
hyp3 <- hypothesis(m3, hyp, class = NULL)

hyp3

plot(density(hyp3$samples$H1), xlim=c(0,1))
abline(v=hyp3$hypothesis$Estimate, col=4, lty=2)
abline(v=hyp3$hypothesis$CI.Lower, col=2, lty=2)
abline(v=hyp3$hypothesis$CI.Upper, col=2, lty=2)


# Looks like lambda is estimated at half of original value with the VCV cor=FALSE

# Why? is this related to the depth of the tree?
V <- vcv(tree, cor=F)
Vcor <- vcv(tree, cor=T)

max(V)
max(Vcor)




# STANDING QUESTIONS

# Does standardizing trait values have same effect?

# Simulate with varying tree depths

# Is phylo hertiability the same as ICC?

# What options do we have extending this to non-Gaussian models?

# Can we use Gaussian kernels to mimic evolutionary models?

