## Create custom Siler distribution

RLW0 <- nimbleRcall(function(x = double(0)){}, Rfun = 'lambertW0',
                   returnType = double(0))

## use optimise function
optimiseNim <- nimbleRcall(function(interval = double(1), 
                                    a1 = double(0),
                                    a2 = double(0), b1 = double(0), b2 = double(0), 
                                    c1 = double(0), p = double(0)){}, Rfun = "optimiseRNim",
                           returnType = double(0))

## R function to run optimiser
optimiseRNim <- function(interval, a1, a2, b1, b2, c1, p) {
  silerRootFn <- function(x, a1, a2, b1, b2, c1, p) {
    abs(pSilerNim(x, a1, a2, b1, b2, c1) - p)
  }
  optimise(silerRootFn, interval, a1, a2, b1, b2, c1, p)$minimum
}


## probability density function
dsilerNim <- nimbleFunction(
  run = function(x = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c1 = double(0), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c1 < 0) {
      return(NaN)
    }
    logS <- (a1 / b1) * (exp(-b1 * x) - 1) - c1 * x + (a2 / b2) * (1 - exp(b2 * x))
    logH <- log(a1 * exp(-b1 * x) + c1 + a2 * exp(b2 * x))
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rsilerNim <- nimbleFunction(
  run = function(n = integer(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c1 = double(0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c1 < 0) {
      return(NaN)
    }
    if(n != 1) print("rsilerNim only allows n = 1; using n = 1.")
    ## sample from two independent distributions and take minimum
    u <- runif(1, 0, 1)
    w0 <- (a1 / c1) * exp((log(u) + a1 / b1) * (b1 / c1))
    w0 <- RLW0(w0)
    x1 <- (-1 / c1) * (log(u) + a1 / b1) + w0 / b1
    x2 <- log(1 - log(runif(1, 0, 1)) * (b2 / a2)) / b2
    xm <- min(x1, x2)
    return(xm)
  })


## cumulative distribution function (and survivor function)
psilerNim <- nimbleFunction(
  run = function(q = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c1 = double(0), 
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c1 < 0) {
      return(NaN)
    }
    logS <- (a1 / b1) * (exp(-b1 * q) - 1) - c1 * q + (a2 / b2) * (1 - exp(b2 * q))
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

## quantile function (not yet working)
qsilerNim <- nimbleFunction(
  run = function(p = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c1 = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c1 < 0) {
      return(NaN)
    }
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    ## find an upper bound
    qU <- 10
    pU <- psilerNim(qU, a1, a2, b1, b2, c1)
    while(qU < p) {
      qU <- qU + 10
      pU <-  psilerNim(qU, a1, a2, b1, b2, c1)
    }
    
    ## run optimiser
    q <- optimiseNim(c(0,qU), a1 = a1, a2 = a2, b1 = b1, b2 = b2, c1 = c1, p = p)
    return(q)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dsilerNim = list(
    BUGSdist = "dsilerNim(a1, a2, b1, b2, c1)",
    Rdist = "dsilerNim(a1, a2, b1, b2, c1)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))
