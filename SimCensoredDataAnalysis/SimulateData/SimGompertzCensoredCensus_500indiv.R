## set up right censored data Gompertz distributed

## set parameters
a <- 0.121
b <- 0.105

nind <- 500

## create sims
sims <- rGompertz(nind, a = a, b = b)

## censor results (1 = interval, 2 = right)
censored <- rbinom(length(sims), 1, censorRate) + 1

## adjust death times to make interval censored
sims <- ceiling(sims)

## adjust sims with random right-censored time
adjust <- which(censored == 2)
sims[adjust] <- sapply(sims[adjust], function(x) {
  sample(1:x, 1)
})

