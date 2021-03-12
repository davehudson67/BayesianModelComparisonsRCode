## load libraries
library(tidyverse)
library(boot)
library(lamW)
library(nimble)
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
logimpweight_s <- readRDS("outputs/logimpweight_s.rds")
logimpweight_ssdA1 <- readRDS("outputs/logimpweight_ssdA1.rds")
logimpweight_ssdA2 <- readRDS("outputs/logimpweight_ssdA2.rds")
logimpweight_ssdB1 <- readRDS("outputs/logimpweight_ssdB1.rds")
logimpweight_ssdB2 <- readRDS("outputs/logimpweight_ssdB2.rds")
logimpweight_ssdC1 <- readRDS("outputs/logimpweight_ssdC1.rds")

## generate log marginal likelihoods
logmarg_s <- log_sum_exp_marg(logimpweight_s)
logmarg_ssdA1 <- log_sum_exp_marg(logimpweight_ssdA1)
logmarg_ssdA2 <- log_sum_exp_marg(logimpweight_ssdA2)
logmarg_ssdB1 <- log_sum_exp_marg(logimpweight_ssdB1)
logmarg_ssdB2 <- log_sum_exp_marg(logimpweight_ssdB2)
logmarg_ssdC1 <- log_sum_exp_marg(logimpweight_ssdC1)

## bootstrap samples
imp_boot_s <- BootsPlot(logimpweight_s, 5000)
imp_boot_ssdA1 <- BootsPlot(logimpweight_ssdA1, 5000)
imp_boot_ssdA2 <- BootsPlot(logimpweight_ssdA2, 5000)
imp_boot_ssdB1 <- BootsPlot(logimpweight_ssdB1, 5000)
imp_boot_ssdB2 <- BootsPlot(logimpweight_ssdB2, 5000)
imp_boot_ssdC1 <- BootsPlot(logimpweight_ssdC1, 5000)

## add prior model weights
priorp <- 1/6
ps <- logmarg_s + log(priorp)
pssdA1 <- logmarg_ssdA1 + log(priorp)
pssdA2 <- logmarg_ssdA2 + log(priorp)
pssdB1 <- logmarg_ssdB1 + log(priorp)
pssdB2 <- logmarg_ssdB2 + log(priorp)
pssdC1 <- logmarg_ssdC1 + log(priorp)

p <- c(ps, pssdA1, pssdA2, pssdB1, pssdB2, pssdC1)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  S = imp_boot_s, 
  sdA1 = imp_boot_ssdA1, 
  sdA2 = imp_boot_ssdA2,
  sdB1 = imp_boot_ssdB1, 
  sdB2 = imp_boot_ssdB2, 
  sdC1 = imp_boot_ssdC1
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
