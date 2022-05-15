# This code was started June 22
# The aim is to load the phylogeny and change the species names to the class names

# increasing the iterations to 4000:3000 - 4 divergent transitions
rm(list = ls()) 
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = FALSE)

library(phytools)
library(ape)
require(rstan)
require(caper)
require(shinystan)
require(reshape2)
library(stringr)
library(ggplot2)
library(plyr)
library(dplyr)

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/pmm/analyses/sync")
} else if(length(grep("Lizzie", getwd()) > 0)) {
  setwd("~/Documents/git/projects/others/deirdre/Synchrony")
} else{
  setwd("/home/deirdre/Synchrony") # for midge
}


##################################################################
# Taking the phylogney code from OSPREE and trying to adapt it for my synchrony project:
#get the synchrony data
#source("Rcode/combiningallphenodata.R")
data <- read.csv("input/synchrony_data.csv")

tree <- read.tree("input/final_tree_pheno.tre")


#get spp and phylo matrix in the same order
head(tree$tip.label)
length(tree$tip.label)

# Full phylogeny
d <- data[match(tree$tip.label, data$phylo.name),]
d$sp.pheno <- str_replace(d$sp.pheno, "__", "_")
d$sp.pheno <- str_replace(d$sp.pheno, " ", "")

#########################################################################################
#### Make a plot for Mammals, birds
class.dat<- read.csv("input/taxonlvl_spname.csv")
class.dat <- class.dat[,c("species.name","class")]

data.temp <- data[,c("species.name","sp.pheno")]
data.temp$count <- 1
data.mini <- aggregate(data.temp["count"], data.temp[c("species.name", "sp.pheno")], FUN = sum)

class <- merge(data.mini, class.dat, by = "species.name", all = T)
sort(unique(class$class)) # "Amphibia" "Aves" "Magnoliopsida" "Mammalia"

amphibia <- subset(class, class == "Amphibia");
amphibia.spp <- unique(amphibia$sp.pheno)

# subset the data and prune the tree:
amphib.dat <- data[data$sp.pheno %in% amphibia.spp,]
amphib_tree <- keep.tip(tree, amphibia.spp)

phymatch.amphib <- data.frame(sp.pheno = amphib_tree$tip.label, sppnum = c(1:length(amphib_tree$tip.label)))
phymatch.amphib$sp.pheno <- str_replace(phymatch.amphib$sp.pheno, "__", "_")
phymatch.amphib$sp.pheno <- str_replace(phymatch.amphib$sp.pheno, " ", "")

d.amphib <- merge(amphib.dat, phymatch.amphib, by="sp.pheno")

d.amphib <- d.amphib[order(d.amphib$sppnum),]
nspecies <- max(d.amphib$sppnum)
#nspecies <- 1275
cophen_tree <- cophenetic(amphib_tree)
vcv_tree <- vcv(amphib_tree, cor = TRUE)

amphib.dat$pheno.fact <- as.numeric(as.factor(amphib.dat$sp.pheno))
amphib.dat$study.fact <- as.numeric(as.factor(amphib.dat$studyid))

length(unique(amphib.dat$sp.pheno)) #1275 w/o thackeray 789

phypriors <- list(
  a_z_prior_mu = 0,
  a_z_prior_sigma = 50,
  astudy_prior_mu = 0,
  astudy_prior_sigma = 50,
  lam_interceptsa_prior_alpha = 1, #
  lam_interceptsa_prior_beta = 1, #
  sigma_interceptsa_prior_mu = 40,
  sigma_interceptsa_prior_sigma = 20,
  b_z_prior_mu = 0,
  b_z_prior_sigma = 10,
  lam_interceptsb_prior_alpha = 2,
  lam_interceptsb_prior_beta = 5, #
  sigma_interceptsb_prior_mu = 5,
  sigma_interceptsb_prior_sigma = 5,
  sigma_y_prior_mu = 20,
  sigma_y_prior_sigma = 10
)


# Function for generating "good" initial values
simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  b_z.temp <- rnorm(n = nspecies, mean = phypriors[["b_z_prior_mu"]], sd = phypriors[["b_z_prior_sigma"]])
  return(append(list(
    a = a_z.temp,
    b_year = b_z.temp),
    phypriors))
}

datalist.amphib <- append(list(N = nrow(amphib.dat),
                        Nspp = nspecies,
                         Nstudy = length(unique(amphib.dat$studyid)),
                        ypred = amphib.dat$doy,
                        species = amphib.dat$pheno.fact,
                         study = amphib.dat$study.fact,
                        year = amphib.dat$yr1980,
                        Vphy = vcv_tree), phypriors)

mdl.amphib <- stan("stan/phylogeny_synchrony_cholesky_grandmean_study.stan",
                 data = datalist.amphib,
                 init = simu_inits,
                 iter = 2000,
                 warmup = 1000,
                 chains = 4,
                 #seed = 62921,
                 refresh = 100
)

save( mdl.amphib, file = "output/phylo_amphibians_grand_study.Rda")

##################################################################
aves <- subset(class, class == "Aves");
aves.spp <- unique(aves$sp.pheno)

aves.dat <- data[data$sp.pheno %in% aves.spp,]
aves_tree <- keep.tip(tree, aves.spp)


phymatch.aves <- data.frame(sp.pheno = aves_tree$tip.label, sppnum = c(1:length(aves_tree$tip.label)))
phymatch.aves$sp.pheno <- str_replace(phymatch.aves$sp.pheno, "__", "_")
phymatch.aves$sp.pheno <- str_replace(phymatch.aves$sp.pheno, " ", "")

d.aves <- merge(aves.dat, phymatch.aves, by="sp.pheno")

d.aves <- d.aves[order(d.aves$sppnum),]
nspecies <- max(d.aves$sppnum)
#nspecies <- 1275
cophen_tree <- cophenetic(aves_tree)
vcv_tree <- vcv(aves_tree, cor = TRUE)

aves.dat$pheno.fact <- as.numeric(as.factor(aves.dat$sp.pheno))
aves.dat$study.fact <- as.numeric(as.factor(aves.dat$studyid))

length(unique(aves.dat$sp.pheno)) #1275 w/o thackeray 789

phypriors <- list(
  a_z_prior_mu = 150,
  a_z_prior_sigma = 50,
  astudy_prior_mu = 0,
  astudy_prior_sigma = 50,
  lam_interceptsa_prior_alpha = 1, #
  lam_interceptsa_prior_beta = 1, #
  sigma_interceptsa_prior_mu = 40,
  sigma_interceptsa_prior_sigma = 20,
  b_z_prior_mu = 0,
  b_z_prior_sigma = 10,
  lam_interceptsb_prior_alpha = 2,
  lam_interceptsb_prior_beta = 5, #
  sigma_interceptsb_prior_mu = 5,
  sigma_interceptsb_prior_sigma = 5,
  sigma_y_prior_mu = 20,
  sigma_y_prior_sigma = 10
)

# Function for generating "good" initial values
simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  b_z.temp <- rnorm(n = nspecies, mean = phypriors[["b_z_prior_mu"]], sd = phypriors[["b_z_prior_sigma"]])
  return(append(list(
    a = a_z.temp,
    b_year = b_z.temp),
    phypriors))
}

datalist.aves <- append(list(N = nrow(aves.dat),
                               Nspp = nspecies,
                               Nstudy = length(unique(aves.dat$studyid)),
                               ypred = aves.dat$doy,
                               species = aves.dat$pheno.fact,
                               study = aves.dat$study.fact,
                               year = aves.dat$yr1980,
                               Vphy = vcv_tree), phypriors)

mdl.aves <- stan("stan/phylogeny_synchrony_cholesky_grandmean_study.stan",
                   data = datalist.aves,
                   init = simu_inits,
                   iter = 4000,
                   warmup = 3000,
                   chains = 4,
                   #seed = 62921,
                   refresh = 20
)

save( mdl.aves, file = "output/phylo_aves_grand_study.Rda")

