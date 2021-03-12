## set up right censored data Exponential distributed

## set parameters
a1 <- exp(-2)
a2 <- exp(-11)
b1 <- 1.2
b2 <- 0.11
c1 <- exp(-5)

## create sims
sims <- rSiler(nind, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c1 = c1)

## censor results (1 = interval, 2 = right)
censored <- rbinom(length(sims), 1, censorRate) + 1

## adjust death times to make interval censored
sims <- ceiling(sims)

## adjust sims with random right-censored time
adjust <- which(censored == 2)
sims[adjust] <- sapply(sims[adjust], function(x) {
  sample(1:x, 1)
})

