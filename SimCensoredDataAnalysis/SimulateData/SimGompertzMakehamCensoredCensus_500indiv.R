## set up right censored data Exponential distributed

## set parameters
a <- 0.0000067
b <- 0.125
c1 <- 0.011

## create sims
sims <- rGompzMake(nind, a = a, b = b, c1 = c1)

inf_error <- which(is.infinite(sims))

if(length(inf_error) > 0){
#while(length(which(is.infinite(sims))) > 0){
  sims[inf_error] <- rGompzMake(length(inf_error), a, b, c1)
}

## adjust sims with random right-censored time
adjust <- which(censored == 2)
sims[adjust] <- sapply(sims[adjust], function(x) {
  sample(1:x, 1)
})

