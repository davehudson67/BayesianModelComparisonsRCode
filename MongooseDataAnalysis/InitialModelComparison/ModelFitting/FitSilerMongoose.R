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
library(survminer)
library(survival)
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
set.seed(seeds[5])

## set up plot output file
pdf("outputs/Siler.pdf")

###########################################################
##                                                      ###
##        Now fit Siler model                           ###
##                                                      ###
###########################################################

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2, c1)
    
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  
})

## set up other components of model
consts <- list(nind = nind)
data <- list(cint = cint, censored = censored, tD = tD)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + runif(1, 0, 1)
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c1 = pars[5], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(5, 100), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    b2 = pars[4],
    c1 = pars[5]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "b2", "c1"), thin = 1)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "b1", "a2", "c1", "b2"), type = 'AF_slice')
#config$addSampler(target = c("a1", "b1"), type = 'AF_slice')

## check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c1"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 5000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## plot mcmcm
plot(run$samples)
samples <- run$samples

## save MCMC
saveRDS(samples, "outputs/samples_s.rds")

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)
preds <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c1 = pars[5], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

## plot ignoring sex differences
pred_plot <- km_plot[[1]] +
  geom_line(aes(x = t, y = Median), data = preds) +
  geom_ribbon(aes(x = t, y = Median, ymin = LCI, ymax = UCI), data = preds, alpha = 0.3)
pred_plot

## pairs plot
samples <- as.matrix(samples)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 10000
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1]
colnames(props) <- c("a1", "a2", "b1", "b2", "c1")

## take random samples from prior (to create defense mixture)
defense <- matrix(rexp(5 * (nimp - nmix), 1), ncol = 5)
colnames(defense) <- c("a1", "a2", "b1", "b2", "c1")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:5)

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights
## log-likelihood function
stopifnot(length(table(censored)) == 2)
log.like <- function(zL, zU, censored, a1, a2, b1, b2, c1) {
  ## calculate log-likelihoods
  llI <- log(pSiler(zU[censored == 1],  a1, a2, b1, b2, c1) - pSiler(zL[censored == 1], a1, a2, b1, b2, c1))
  llR <- pSiler(zL[censored == 2], a1, a2, b1, b2, c1, log = TRUE, lower.tail = FALSE)
  sum(llI) + sum(llR)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, zL, zU, censored) {
                           log.like(zL, zU, censored, pars[1], pars[2], pars[3], pars[4], pars[5])
                         }, zL = zL, zU = zU, censored = censored, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) +
  dexp(props[, 4], 1, log = TRUE) + dexp(props[, 5], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) + dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) + dexp(props[, 4], 1, log = TRUE) +
                                                                              dexp(props[, 5], 1, log = TRUE)))
saveRDS(logimpweight, "outputs/logimpweight_s.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
BootsPlot(logimpweight, 5000, trace = TRUE)

## turn graphics device off
dev.off()
