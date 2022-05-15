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

# Amphibians:
amphibia <- subset(phylo_pheno, class == "Amphibia");

amphibia.fam <- c("Hylidae", "Ambystomatidae", "Bufonidae", "Plethodontidae", "Microhylidae", "Hynobiidae", "Ranidae", "Rhacophoridae", "Pelobatidae", "Salamandridae")
amphib_tree <- keep.tip(tree, amphibia.fam)

d <- phylo_pheno[match(amphib_tree$tip.label, phylo_pheno$family),]; length(unique(d$family))

phymatch <- data.frame(family = amphib_tree$tip.label, sppnum = c(1:length(amphib_tree$tip.label)))

d <- merge(phylo_pheno, phymatch, by="family")

length(unique(d$sp.pheno)) # 34
length(unique(d$studyid)) # 14
length(unique(d$species.name)) # 34

spt <- unique(d$sp.pheno)

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)

cophen_tree <- cophenetic(amphib_tree)
vcv_tree <- vcv(amphib_tree, cor = TRUE)

amphibia$pheno.fact <- as.numeric(as.factor(as.character(amphibia$sp.pheno)))
amphibia$species.fact <- as.numeric(as.factor(as.character(amphibia$species.name)))
amphibia$family.fact <- as.numeric(as.factor(as.character(amphibia$family)))
amphibia$study.fact <- as.numeric(as.factor(as.character(amphibia$studyid)))


phypriors <- list(
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
   a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  b_z.temp <- rnorm(n = nspecies, mean = phypriors[["b_z_prior_mu"]], sd = phypriors[["b_z_prior_sigma"]])
  return(append(list(
    b_year = b_z.temp),
    phypriors))
}

Jtemp <- amphibia[,c("pheno.fact","family.fact")]
Jtemp$count <- 1
J <-  aggregate(Jtemp["count"], Jtemp[c("pheno.fact","family.fact")], FUN = sum)

datalist <- append(list(N = nrow(amphibia),
                        Nspp = length(unique(amphibia$sp.pheno)),
                        Nfam = length(unique(amphibia$family)),
                        sppnum = amphibia$pheno.fact,
                        famnum = J$family.fact,
                        y = amphibia$doy,
                        x = amphibia$yr1980,
                        Vphy = vcv_tree), phypriors)

sort(amphibia$pheno.fact)
mdl.fam.amphib <- stan("..//stan/phylogeny_family_cholesky.stan",
                data = datalist,
               # init = simu_inits,
                iter = 4000,
                warmup = 2000,
                chains = 4,
                seed = 62921,
                refresh = 10
)


save( mdl.fam.amphib, file = "output/fam_phylo_nostudy.Rda")

load("output/fam_phylo_nostudy.Rda")

ssm <-  as.shinystan(mdl.fam.amphib)
launch_shinystan(ssm)
##########################################################

# Aves

aves <- subset(phylo_pheno, class == "Aves")
aves.fam <- as.character(unique(aves$family))

aves_tree <- keep.tip(tree, aves.fam)

d <- aves[match(aves_tree$tip.label, aves$family),]; length(unique(aves$family))

phymatch <- data.frame(family = aves_tree$tip.label, sppnum = c(1:length(aves_tree$tip.label)))

d <- merge(aves, phymatch, by="family")

spt <- unique(d$sp.pheno)

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 1275
cophen_tree <- cophenetic(aves_tree)
vcv_tree <- vcv(aves_tree, cor = TRUE)

aves$pheno.fact <- as.numeric(as.factor(as.character(aves$sp.pheno)))
aves$species.fact <- as.numeric(as.factor(as.character(aves$species.name)))
aves$family.fact <- as.numeric(as.factor(as.character(aves$family)))
aves$study.fact <- as.numeric(as.factor(as.character(aves$studyid)))


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

#phylo_pheno <- phylo_pheno[complete.cases(phylo_pheno$doy),]


Jtemp <- aves[,c("pheno.fact","family.fact")]
Jtemp$count <- 1
J <-  aggregate(Jtemp["count"], Jtemp[c("pheno.fact","family.fact")], FUN = sum)

datalist <- append(list(N = nrow(aves),
                        Nspp = length(unique(aves$sp.pheno)),
                        Nfam = length(unique(aves$family)),
                        sppnum = aves$pheno.fact,
                        famnum = J$family.fact,
                        y = aves$doy,
                        # study = phylo_pheno$study.fact,
                        x = aves$yr1980,
                        Vphy = vcv_tree), phypriors)
datalist$famnum


mdl.fam <- stan("..//stan/phylogeny_family_cholesky.stan",
                data = datalist,
                #init = simu_inits,
                iter = 4000,
                warmup = 2000,
                chains = 1,
                seed = 62921,
                refresh = 10
)


