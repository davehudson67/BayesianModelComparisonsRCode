## set up right censored data Exponential distributed

## set parameters
r <- 0.005
censorRate <- 0.1
nind <- 500

## create sims
sims <- rexp(nind, r = r)

## censor results (1 = interval, 2 = right)
censored <- rbinom(length(sims), 1, censorRate) + 1

## adjust death times to make interval censored
sims <- ceiling(sims)

## adjust sims with random right-censored time
adjust <- which(censored == 2)
sims[adjust] <- sapply(sims[adjust], function(x) {
  sample(1:x, 1)
})

