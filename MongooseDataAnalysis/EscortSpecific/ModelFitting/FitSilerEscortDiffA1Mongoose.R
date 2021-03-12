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
set.seed(seeds[3])

## set up plot output file
pdf("outputs/SilerEscortDiffA1.pdf")

###########################################################
##                                                      ###
##   Now fit Siler model with escort differences on a1  ###
##                                                      ###
###########################################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1[escort[i] + 1], a2, b1, b2, c1)
  }
  
  for (j in 1:g){
    ## priors
    a1[j] ~ dexp(1)
  }
  
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)

})

## set up other components of model
consts <- list(nind = nind, g = g, escort = escort)

data <- list(cint = cint, censored = censored, tD = tD)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + 1
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored, escort) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, escort) {
    if(any(pars < 0)) {
      return(NA)
    }
    llL <- sum(dSiler(t[escort == 0], a1 = pars[1], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
    llH <- sum(dSiler(t[escort == 1], a1 = pars[2], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
    llL + llH
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    ## optimise to interval-censored only
    pars <- optim(rexp(6, 100), optFn, t = tD, escort = escort, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  ## check log-likelihoods
  ll <- sum(dSiler(tD[escort == 0], a1 = pars[1], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
  ll <- ll + sum(dSiler(tD[escort == 1], a1 = pars[2], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
  stopifnot(is.finite(ll))
  ## output initial values
  list(
    tD = tD,
    a1 = c(pars[1], pars[2]),
    a2 = pars[3], 
    b1 = pars[4],
    b2 = pars[5],
    c1 = pars[6]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored, escort))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("a1", "a2", "b1", "b2", "c1"), thin = 1)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "a2", "b1", "b2", "c1"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c1"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 2000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary

#Plot mcmcm
samples <- run$samples
plot(run$samples)

## save MCMC
saveRDS(samples, "outputs/post_sedA1.rds")

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)

predsL <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

predsH <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[2], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

## plot posterior predictive K-M plots against the posterior predictive survival curve
pred_escort <- mclapply(as.matrix(samples)[1:2000, 1], function(p, temp_dat) {
  temp_dat$escort <- factor(as.character(temp_dat$escort_cat))
  temp_fit <- surv_fit(Surv(tD, cens) ~ escort, data = temp_dat)
  tibble(time = temp_fit$time, surv = temp_fit$surv, escort = rep(names(temp_fit$strata), temp_fit$strata))
}, temp_dat = mong_km_dat, mc.cores = 20) %>%
  bind_rows() %>%
  group_by(time, escort) %>%
  summarise(
    Median = median(surv),
    LCI = quantile(surv, probs = 0.025),
    UCI = quantile(surv, probs = 0.975)
  )
pred_plot <- ggplot(pred_escort, aes(x = time)) +
  geom_point(aes(y = Median, colour = escort)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = escort), alpha = 0.3) +
  geom_line(aes(x = t, y = Median), data = predsL, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = predsL, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = predsL, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = Median), data = predsH, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = predsH, linetype = "dotted", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = predsH, linetype = "dotted", alpha = 0.7) +
  xlab("Time") + ylab("Survival") + labs(colour = "Escort", fill = "Escort") + 
  ggtitle("Comparison Survival Curves", subtitle = "KM plots of actual data against predicted escort specific differences in a1")
pred_plot

## pairs plot
samples <- as.matrix(samples)
colnames(samples) <- c("a1l", "a1h", "a2", "b1", "b2", "c1")
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
nimp <- 500000
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1]
colnames(props) <- c("a1l", "a1h", "a2", "b1", "b2", "c1")

## take random samples from prior (to create defense mixture)
defense <- matrix(rexp(6 * (nimp - nmix), 1), ncol = 6)
colnames(defense) <- c("a1l", "a1h", "a2", "b1", "b2", "c1")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:6)

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights

## log-likelihood function
## REQUIRES BOTH RIGHT- AND INTERVAL-CENSORED INDIVIDUALS
## AND MISSING SEX INDIVIDUALS
log.like <- function(zL, zU, censored, escort, pars) {
  ## calculate log-likelihoods
  ## interval-censored
  llI_l <- log(pSiler(zU[censored == 1 & escort == 0], pars[1], pars[3], pars[4], pars[5], pars[6]) - pSiler(zL[censored == 1 & escort == 0], pars[1], pars[3], pars[4], pars[5], pars[6]))
  llI_h <- log(pSiler(zU[censored == 1 & escort == 1], pars[2], pars[3], pars[4], pars[5], pars[6]) - pSiler(zL[censored == 1 & escort == 1], pars[2], pars[3], pars[4], pars[5], pars[6]))
  
  ## right-censored
  llR_l <- pSiler(zL[censored == 2 & escort == 0], pars[1], pars[3], pars[4], pars[5], pars[6], log = TRUE, lower.tail = FALSE)
  llR_h <- pSiler(zL[censored == 2 & escort == 1], pars[2], pars[3], pars[4], pars[5], pars[6], log = TRUE, lower.tail = FALSE)
  
    ## return log-likelihood
  sum(llI_l) + sum(llR_h) + sum(llI_l) + sum(llR_h)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, zL, zU, censored, escort) {
                           log.like(zL, zU, censored, escort, pars)
                         }, zL = zL, zU = zU, censored = censored, escort = escort, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) +
  dexp(props[, 4], 1, log = TRUE) + dexp(props[, 5], 1, log = TRUE) +
  dexp(props[, 6], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) + dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) + dexp(props[, 4], 1, log = TRUE) +
                                                                              dexp(props[, 5], 1, log = TRUE) + dexp(props[, 6], 1, log = TRUE)))
saveRDS(logimpweight, "outputs/logimpweight_sedA1.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
BootsPlot(logimpweight, 5000, trace = TRUE)

## turn graphics device off
dev.off()
