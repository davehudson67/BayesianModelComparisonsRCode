###################################################
##                                               ##
##   refit the data using an exponential model   ##
##                                               ##
###################################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
  ## likelihood for interval-truncated exponential
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dexp(r)
  }
  
  ## priors
  r ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
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
initFn <- function(cint, censored) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, censored) {
    if(any(pars < 0)) {
      return(NA)
    }
    ## fit to just interval-censored
    ll <- sum(dexp(t[censored == 1], r = pars[1], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    ## optimise to interval-censored only
    pars <- optim(rexp(10, 10), optFn, t = tD, censored = censored, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  ## set suitable initial values for right-censored
  tD[censored == 2] <- rExpo(sum(censored == 2), r = pars[1], zL = cint[censored == 2, 2])
  
  ## output initial values
  list(
    tD = tD,
    r = pars[1]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("r"), thin = 1)

## check monitors and samplers
config$printMonitors()
config$printSamplers("r")

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

## pairs plot
samples <- as.matrix(samples)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), , drop = FALSE]
samples %>%
  as.data.frame() %>%
  ggplot(aes(x = r)) + geom_density()

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 10000
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- "r"

## take random samples from prior (to create defense mixture)
defense <- matrix(rexp(1 * (nimp - nmix), 1), ncol = 1)

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggplot(aes(x = r, colour = type)) + geom_density()

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights
## log-likelihood function
log.like <- function(zL, zU, censored, r) {
  ## calculate log-likelihoods
  llI <- log(pexp(zU[censored == 1], r) - pexp(zL[censored == 1], r))
  llR <- pexp(zL[censored == 2], r, log = TRUE, lower.tail = FALSE)
  sum(llI) + sum(llR)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
    function(pars, zL, zU, censored) {
        log.like(zL, zU, censored, pars)
    }, zL = zL, zU = zU, censored = censored, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE)))
saveRDS(logimpweight, "logimpweight_e.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## log marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

#Bootstrap the importance weights to create 95% intervals
imp_boot <- BootsPlot(logimpweight, 5000)
