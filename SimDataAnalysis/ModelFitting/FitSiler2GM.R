###########################################################
##                                                      ###
##        Now fit Siler model                           ###
##                                                      ###
###########################################################

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated Siler
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
data <- list(tD = tD)
initFn <- function(tD) {
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
      pars <- optim(rexp(5, 100), optFn, t = tD, control = list(fnscale = -1))
      k <- k + 1
    }
    if(k == 20) {
      stop("Can't sample initial values")
    }
    pars <- pars$par
    list(
      a1 = pars[1], 
      a2 = pars[2], 
      b1 = pars[3],
      b2 = pars[4],
      c1 = pars[5]
    )
  }
  
inits <- initFn(tD)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(tD))

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "b2", "c1"), thin = 1)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "a2","b1", "b2", "c1"), type = 'AF_slice')
#config$addSampler(target = c(), type = 'AF_slice')

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
log.like <- function(x, a1, a2, b1, b2, c1) {
  sum(dSiler(x, a1, a2, b1, b2, c1, log = TRUE))
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, x) {
                           log.like(x, pars[1], pars[2], pars[3], pars[4], pars[5])
                         }, x = tD, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) +
  dexp(props[, 4], 1, log = TRUE) + dexp(props[, 5], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) +
                                                                              dexp(props[, 4], 1, log = TRUE) + dexp(props[, 5], 1, log = TRUE)))
saveRDS(logimpweight, "logimpweight_s.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
imp_boot <- BootsPlot(logimpweight, 5000)
