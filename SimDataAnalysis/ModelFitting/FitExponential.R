##################################
##                              ## 
##  Now fit Exponential model   ##
##                              ##
##################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
  ## likelihood for interval-truncated Siler
    tD[i] ~ dexp(r)
  }
  
  ## priors
  r ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(tD = tD)
inits <- list(
  r = rexp(1, 100) 
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = "r", thin = 1)

## check monitors and samplers
config$printMonitors()
config$printSamplers("r")

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 10000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## plot mcmcm
plot(run$samples)
samples <- run$samples

## 
samples <- as.matrix(samples)
colnames(samples) <- "r"

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
colnames(props) <- c("r")

## take random samples from prior (to create defense mixture)
defense <- matrix(rexp(1 * (nimp - nmix), 1), ncol = 1)
colnames(defense) <- c("r")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:1)

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights
## log-likelihood function
log.like <- function(x, r) {
  sum(dexp(x, r, log = TRUE))
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, x) {
                           log.like(x, pars[1])
                         }, x = tD, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE)))
saveRDS(logimpweight, "logimpweight_e.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
imp_boot <- BootsPlot(logimpweight, 5000)
