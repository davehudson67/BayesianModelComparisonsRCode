## load libraries
library(tidyverse)
library(boot)
library(lamW)
rm(list = ls())

## source necessary R distributions
source("../../Distributions/Dist_Gompertz.R")
source("../../Distributions/Dist_GompertzNim.R")
source("../../Distributions/Dist_GompertzMakeham.R")
source("../../Distributions/Dist_GompertzMakehamNim.R")
source("../../Distributions/Dist_Siler.R")
source("../../Distributions/Dist_SilerNim.R")
source("../../Distributions/Dist_Expo.R")

## source additional R functions
source("../../ModelComparison_FUNCTIONS.R")

## load data
load("mong.RData")

## set seed according to model
set.seed(seeds[7])

## set up plot output file
pdf("outputs/ModelComparisons.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples
logimpweight_sned <- readRDS("outputs/logimpweight_sned.rds")
logimpweight_sedA1 <- readRDS("outputs/logimpweight_sedA1.rds")
logimpweight_sedA2 <- readRDS("outputs/logimpweight_sedA2.rds")
logimpweight_sedB1 <- readRDS("outputs/logimpweight_sedB1.rds")
logimpweight_sedB2 <- readRDS("outputs/logimpweight_sedB2.rds")
logimpweight_sedC1 <- readRDS("outputs/logimpweight_sedC1.rds")

## generate log marginal likelihoods
logmarg_sned <- log_sum_exp_marg(logimpweight_sned)
logmarg_sedA1 <- log_sum_exp_marg(logimpweight_sedA1)
logmarg_sedA2 <- log_sum_exp_marg(logimpweight_sedA2)
logmarg_sedB1 <- log_sum_exp_marg(logimpweight_sedB1)
logmarg_sedB2 <- log_sum_exp_marg(logimpweight_sedB2)
logmarg_sedC1 <- log_sum_exp_marg(logimpweight_sedC1)

## bootstrap samples
imp_boot_sned <- BootsPlot(logimpweight_sned, 5000)
imp_boot_sedA1 <- BootsPlot(logimpweight_sedA1, 5000)
imp_boot_sedA2 <- BootsPlot(logimpweight_sedA2, 5000)
imp_boot_sedB1 <- BootsPlot(logimpweight_sedB1, 5000)
imp_boot_sedB2 <- BootsPlot(logimpweight_sedB2, 5000)
imp_boot_sedC1 <- BootsPlot(logimpweight_sedC1, 5000)

## add prior model weights
priorp <- 1/6
psned <- logmarg_sned + log(priorp)
psedA1 <- logmarg_sedA1 + log(priorp)
psedA2 <- logmarg_sedA2 + log(priorp)
psedB1 <- logmarg_sedB1 + log(priorp)
psedB2 <- logmarg_sedB2 + log(priorp)
psedC1 <- logmarg_sedC1 + log(priorp)

p <- c(psned, psedA1, psedA2, psedB1, psedB2, psedC1)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  S = imp_boot_sned, 
  sedA1 = imp_boot_sedA1, 
  sedA2 = imp_boot_sedA2,
  sedB1 = imp_boot_sedB1, 
  sedB2 = imp_boot_sedB2, 
  sedC1 = imp_boot_sedC1
)
MargLike.plot(mods)

## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg

## turn graphics device off
dev.off()
