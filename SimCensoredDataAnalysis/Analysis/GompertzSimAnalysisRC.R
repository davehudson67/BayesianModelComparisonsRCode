## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(GGally)
library(coda)
library(mclust)
library(parallel)

## set seed
set.seed(41)

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

## simulate right censored data
source("../SimulateData/SimGompertzCensoredCensus_500indiv.R")

## setup data
tD <- rep(NA, length(sims))
tL <- sims - 1
tU <- sims

## define censoring matrices
## (here tL is lower bound and tU upper bound
## and sims is censoring time for right-censored
## so no further adjustments needed)
cint <- cbind(tL, tU)
colnames(cint) <- NULL
cint[censored == 2, 1] <- 0

## set up zL/zU input matrix
zL <- tL
zL[censored == 2] <- tU[censored == 2]

zU <- tU
zU[censored == 2] <- NA

## set maximum observation time
tM <- rep(ceiling(max(sims) + 1), length(sims))

## fit Siler model
source("../ModelFitting/FitSilerCensoredDEFENSE.R")

## fit Gompertz-Makeham model
source("../ModelFitting/FitGompertzMakehamCensoredDEFENSE.R")

## fit Gompertz model
source("../ModelFitting/FitGompertzCensoredDEFENSE.R")

## fit Exponential model
source("../ModelFitting/FitExponentialCensoredDEFENSE.R")

## load IS samples
logimpweight_e <- readRDS("logimpweight_e.rds")
logimpweight_g <- readRDS("logimpweight_g.rds")
logimpweight_gm <- readRDS("logimpweight_gm.rds")
logimpweight_s <- readRDS("logimpweight_s.rds")

## generate log marginal likelihoods
logmarg_e <- log_sum_exp_marg(logimpweight_e)
logmarg_g <- log_sum_exp_marg(logimpweight_g)
logmarg_gm <- log_sum_exp_marg(logimpweight_gm)
logmarg_s <- log_sum_exp_marg(logimpweight_s)

## bootstrap samples
imp_boot_e <- BootsPlot(logimpweight_e, 5000, TRUE)
imp_boot_g <- BootsPlot(logimpweight_g, 5000, TRUE)
imp_boot_gm <- BootsPlot(logimpweight_gm, 5000, TRUE)
imp_boot_s <- BootsPlot(logimpweight_s, 5000, TRUE)

## add prior model weights
pe <- logmarg_e + log(1/4)
pgm <- logmarg_gm + log(1/4)
pg <- logmarg_g + log(1/4)
ps <- logmarg_s + log(1/4)
p <- c(pe, pgm, pg, ps)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  E = imp_boot_e, 
  GM = imp_boot_gm, 
  G = imp_boot_g,
  S = imp_boot_s 
)
MargLike.plot(mods)

## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg
