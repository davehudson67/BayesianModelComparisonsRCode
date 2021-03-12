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
set.seed(seeds[4])

## set up plot output file
pdf("outputs/SilerSexSpecificA2.pdf")

###########################################################
##                                                      ###
##   Now fit Siler model with sex differences on a2     ###
##                                                      ###
###########################################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2[sex[i] + 1], b1, b2, c1)
    sex[i] ~ dbern(Pm)
  }
  
  for (j in 1:g){
    ## priors
    a2[j] ~ dexp(1)
  }
  
  a1 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  Pm ~ dunif(0,1)
  
})

## set up other components of model
consts <- list(nind = nind, g = g)

data <- list(cint = cint, censored = censored, tD = tD, sex = sex)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + rexp(1, 1)
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored, sex) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, sex) {
    if(any(pars < 0)) {
      return(NA)
    }
    llM <- sum(dSiler(t[sex == 0], a1 = pars[1], a2 = pars[2], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
    llF <- sum(dSiler(t[sex == 1], a1 = pars[1], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
    llM + llF
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    ## sample sex proportion
    Pm <- runif(1, 0.4, 0.6)
    ## sample missing sex indicators
    sexI <- rbinom(length(censored), size = 1, prob = Pm)
    sexI[!is.na(sex)] <- sex[!is.na(sex)]
    ## optimise to interval-censored only
    pars <- optim(rexp(6, 100), optFn, t = tD, sex = sexI, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  ## check log-likelihoods
  ll <- sum(dSiler(tD[sexI == 0], a1 = pars[1], a2 = pars[2], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
  ll <- ll + sum(dSiler(tD[sexI == 1], a1 = pars[1], a2 = pars[3], b1 = pars[4], b2 = pars[5], c1 = pars[6], log = TRUE))
  stopifnot(is.finite(ll))
  ## reformat sex initial conditions correctly
  sexI[!is.na(sex)] <- NA
  ## output initial values
  list(
    tD = tD,
    sex = sexI,
    Pm = Pm,
    a1 = pars[1],
    a2 = c(pars[2], pars[3]), 
    b1 = pars[4],
    b2 = pars[5],
    c1 = pars[6]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored, sex))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("a1", "a2", "b1", "b2", "c1", "Pm"), thin = 1)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "b1"), type = 'AF_slice')
config$addSampler(target = c("a2", "c1"), type = 'AF_slice')
#config$addSampler(target = c("a2", "b2", "c1"), type = 'AF_slice')
config$addSampler(target = c("b2"), type = 'slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c1", "Pm"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
    niter = 50000, 
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
saveRDS(samples, "outputs/post_ssda2.rds")

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)

predsM <- as.matrix(samples)[1:2000, -1] %>%
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

predsF <- as.matrix(samples)[1:2000, -1] %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[4], b2 = pars[5], c1 = pars[6], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

## plot posterior predictive K-M plots against the posterior predictive survival curve
pred_sex <- mclapply(as.matrix(samples)[1:2000, 1], function(p, temp_dat) {
  temp_dat$sex[temp_dat$sex == "P"] <- ifelse(rbinom(sum(temp_dat$sex == "P"), size = 1, prob = p) == 1, "M", "F")
  temp_dat$sex <- factor(as.character(temp_dat$sex))
  temp_fit <- surv_fit(Surv(tD, cens) ~ sex, data = temp_dat)
  tibble(time = temp_fit$time, surv = temp_fit$surv, sex = rep(names(temp_fit$strata), temp_fit$strata))
}, temp_dat = mong_km_dat, mc.cores = 20) %>%
  bind_rows() %>%
  mutate(sex = ifelse(sex == "sex=F", "F", "M")) %>%
  group_by(time, sex) %>%
  summarise(
    Median = median(surv),
    LCI = quantile(surv, probs = 0.025),
    UCI = quantile(surv, probs = 0.975)
  )
pred_plot <- ggplot(pred_sex, aes(x = time)) +
  geom_point(aes(y = Median, colour = sex)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = sex), alpha = 0.3) +
  geom_line(aes(x = t, y = Median), data = predsM, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = predsM, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = predsM, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = Median), data = predsF, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = predsF, linetype = "dotted", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = predsF, linetype = "dotted", alpha = 0.7) +
  xlab("Time") + ylab("Survival") + labs(colour = "Sex", fill = "Sex") + 
  ggtitle("Comparison Survival Curves", subtitle = "KM plots of actual data against predicted sex specific differences in a2")
pred_plot

## pairs plot
samples <- as.matrix(samples)
colnames(samples) <- c("Pm", "a1", "a2f", "a2m", "b1", "b2", "c1")
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
colnames(props) <- c("Pm", "a1", "a2f", "a2m", "b1", "b2", "c1")

## take random samples from prior (to create defense mixture)
defense <- matrix(runif(nimp - nmix, 0, 1), ncol = 1)
defense <- cbind(defense, matrix(rexp(6 * (nimp - nmix), 1), ncol = 6))
colnames(defense) <- c("Pm", "a1", "a2f", "a2m", "b1", "b2", "c1")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:7)

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights

## log-likelihood function
## REQUIRES BOTH RIGHT- AND INTERVAL-CENSORED INDIVIDUALS
## AND MISSING SEX INDIVIDUALS
stopifnot(all(dim(table(sex, censored)) == 2))
stopifnot(length(table(is.na(sex))) == 2)
log.like <- function(zL, zU, censored, sex, pars) {
  ## calculate log-likelihoods
  
  ## interval-censored, known sex
  llI_f <- log(pSiler(zU[censored == 1 & sex == 0 & !is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7]) - pSiler(zL[censored == 1 & sex == 0 & !is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7])) + log(1 - pars[1])
  llI_m <- log(pSiler(zU[censored == 1 & sex == 1 & !is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7]) - pSiler(zL[censored == 1 & sex == 1 & !is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7])) + log(pars[1])
  
  ## right-censored, known sex
  llR_f <- pSiler(zL[censored == 2 & sex == 0 & !is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7], log = TRUE, lower.tail = FALSE) + log(1 - pars[1])
  llR_m <- pSiler(zL[censored == 2 & sex == 1 & !is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7], log = TRUE, lower.tail = FALSE) + log(pars[1])
  
  ## interval-censored unknown sex
  llI_miss <- (1 - pars[1]) * (pSiler(zU[censored == 1 & is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7]) - pSiler(zL[censored == 1 & is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7]))
  llI_miss <- llI_miss + pars[1] * (pSiler(zU[censored == 1 & is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7]) - pSiler(zL[censored == 1 & is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7]))
  llI_miss <- log(llI_miss)
  
  ## right-censored unknown sex
  llR_miss <- (1 - pars[1]) * pSiler(zL[censored == 2 & is.na(sex)], pars[2], pars[3], pars[5], pars[6], pars[7], lower.tail = FALSE)
  llR_miss <- llR_miss + pars[1] * pSiler(zL[censored == 2 & is.na(sex)], pars[2], pars[4], pars[5], pars[6], pars[7], lower.tail = FALSE)
  llR_miss <- log(llR_miss)
  
  ## return log-likelihood
  sum(llI_f) + sum(llR_f) + sum(llI_m) + sum(llR_m) + sum(llI_miss) + sum(llR_miss)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
    function(pars, zL, zU, censored, sex) {
        log.like(zL, zU, censored, sex, pars)
    }, zL = zL, zU = zU, censored = censored, sex = sex, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dunif(props[, 1], 0, 1, log = TRUE)
for(j in 2:ncol(props)) {
  logimpweight <- logimpweight + dexp(props[, j], 1, log = TRUE)
}

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dunif(props[, 1], 0, 1, log = TRUE) +
  dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) + 
  dexp(props[, 4], 1, log = TRUE) + dexp(props[, 5], 1, log = TRUE) +
  dexp(props[, 6], 1, log = TRUE) + dexp(props[, 7], 1, log = TRUE)))
saveRDS(logimpweight, "outputs/logimpweight_ssdA2.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
BootsPlot(logimpweight, 5000, trace = TRUE)

## turn graphics device off
dev.off()