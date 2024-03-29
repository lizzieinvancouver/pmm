# code started Nov 23, 2021 by DL
# aim of this code is to see if the model will run on a family level, linnean tree
# major issue I am currenly having is that the for loop the calculates the covariance matrix is taking too long to run on very large matrcies 

# modified January 19, 2022
# The model for the species level is still having issues with it's run time, so we are moving on with a family level tree

rm(list = ls()) 
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = FALSE)

library(phytools)
library(ape)
require(rstan)
require(shinystan)
require(reshape2)
library(stringr)
require(tidybayes)

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/pmm/analyses/sync/familyLevel")
} else if(length(grep("Lizzie", getwd()) > 0)) {
  setwd("~/Documents/git/projects/others/deirdre/Synchrony")
} else{
  setwd("/home/deirdre/pmm") # for midge
}

dat.fin <- read.csv("..//input/synchrony_data.csv")

# ################################################
phylo <- read.csv("..//input/taxonlvl_spname.csv")

phylo$species.name <- as.factor(phylo$species.name)
phylo$species.y <- as.factor(phylo$sp)
phylo$genus <- as.factor(phylo$genus)
phylo$family <- as.factor(phylo$family)
phylo$order <- as.factor(phylo$order)
phylo$class <- as.factor(phylo$class)
phylo$phylum <- as.factor(phylo$phylum)
phylo$kingdom <- as.factor(phylo$kingdom)

phylo_pheno <- merge(phylo, dat.fin, by ="species.name")
# ################################################
tree <- read.tree("..//input/family_lvl_tree.tre")
head(tree$tip.label)
length(tree$tip.label) #217

################################################

d <- phylo_pheno[match(tree$tip.label, phylo_pheno$family),]; length(unique(d$family))

phymatch <- data.frame(family = tree$tip.label, sppnum = c(1:length(tree$tip.label)))

d <- merge(phylo_pheno, phymatch, by="family")

length(unique(d$sp.pheno)) # 1279
length(unique(d$studyid)) # 147
length(unique(d$species.name)) # 1200

spt <- unique(d$sp.pheno)

temp <- phylo_pheno[!phylo_pheno$sp.pheno %in% spt, ]
d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 1275
cophen_tree <- cophenetic(tree)
vcv_tree <- vcv(tree, cor = TRUE)

phylo_pheno$pheno.fact <- as.numeric(as.factor(phylo_pheno$sp.pheno))
phylo_pheno$species.fact <- as.numeric(as.factor(phylo_pheno$species.name))
phylo_pheno$family.fact <- as.numeric(as.factor(phylo_pheno$family))
phylo_pheno$study.fact <- as.numeric(as.factor(phylo_pheno$studyid))


phypriors <- list(
  #a_z_prior_mu = 150,
                  #a_z_prior_sigma = 50,
                  #astudy_prior_mu = 0,
                  #astudy_prior_sigma = 50,
                  # lam_interceptsa_prior_alpha = 1, # 
                  # lam_interceptsa_prior_beta = 1, # 
                  # sigma_interceptsa_prior_mu = 40,
                  # sigma_interceptsa_prior_sigma = 20,
                  b_z_prior_mu = 0,
                  b_z_prior_sigma = 10,
                  lam_interceptsb_prior_alpha = 1, #
                  lam_interceptsb_prior_beta = 1, # 
                  sigma_interceptsb_prior_mu = 5,
                  sigma_interceptsb_prior_sigma = 5,
                  sigma_y_mu_prior = 20,
                  sigma_y_sigma_prior = 10
)

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
  # a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  b_z.temp <- rnorm(n = nspecies, mean = phypriors[["b_z_prior_mu"]], sd = phypriors[["b_z_prior_sigma"]])
  return(append(list(
                     b_year = b_z.temp),
                phypriors))
}

phylo_pheno <- phylo_pheno[complete.cases(phylo_pheno$doy),]

length(unique(phylo_pheno$species.name))
length(unique(phylo_pheno$studyid))
length(unique(phylo_pheno$sp.pheno))
length(unique(phylo_pheno$family))

Jtemp <- phylo_pheno[,c("pheno.fact","family.fact")]
Jtemp$count <- 1
J <-  aggregate(Jtemp["count"], Jtemp[c("pheno.fact","family.fact")], FUN = sum)

datalist <- append(list(N = nrow(phylo_pheno),
                        Nspp = length(unique(phylo_pheno$sp.pheno)),
                        Nfam = length(unique(phylo_pheno$family)),
                        sppnum = phylo_pheno$pheno.fact,
                        famnum = J$family.fact,
                        y = phylo_pheno$doy,
                       # study = phylo_pheno$study.fact,
                        x = phylo_pheno$yr1980,
                        Vphy = vcv_tree), phypriors)


mdl.fam <- stan("..//stan/phylogeny_family_cholesky.stan",
                 data = datalist,
                 init = simu_inits,
                 iter = 4000,
                 warmup = 2000,
                 chains = 4,
                 seed = 62921,
                 refresh = 10
)


save( mdl.fam, file = "output/fam_phylo_nostudy.Rda")


#load("output/fam_phylo_nostudy.Rda")
# # 
ssm <-  as.shinystan(mdl.fam)
launch_shinystan(ssm)
# #
get_variables(mdl.fam)

summary(mdl.fam)$summary[c("b_z","lam_interceptsb", "sigma_interceptsb","sig_a","sig_b","sig_y","sig_a_fam","sig_b_fam","mu_a","mu_b"),"mean"]

pairs(mdl.fam, pars = c("b_z","lam_interceptsb", "sigma_interceptsb"))
