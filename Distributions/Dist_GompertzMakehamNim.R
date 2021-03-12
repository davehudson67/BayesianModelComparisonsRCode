#Create custom Gompertz Makeham distribution:

RLW0<- nimbleRcall(function(x = double(0)){}, Rfun = 'lambertW0',
                   returnType = double(0))

## probability density function
dgompzMakeNim <- nimbleFunction(
  run = function(x = double(0), a = double(0),
                 b = double(0), c1 = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0 | c1 < 0) {
      return(NaN)
    }
    logS <- -c1 * x - (a / b) * (exp(b * x) - 1)
    logH <- log(a * exp(b * x) + c1)
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rgompzMakeNim <- nimbleFunction(
  run = function(n = integer(0), a = double(0),
                 b = double(0), c1 = double(0)) {
    returnType(double(0))
    if(a < 0 | b < 0 | c1 < 0) {
      return(NaN)
    }
    if(n != 1) print("rgompzMakeNim only allows n = 1; using n = 1.")
    u <- runif(1, 0, 1)
    w0 <- (a / c1) * exp(a / c1) * ((1 - u)^(-b / c1))
    w0 <- RLW0(w0)
    rs <- (a / (b * c1)) - (1 / c1) * log(1 - u) - (1 / b) *  w0
    return(rs) 
  })

## cumulative distribution function (and survivor function)
pgompzMakeNim <- nimbleFunction(
  run = function(q = double(0), a = double(0),
                 b = double(0), c1 = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0 | c1 < 0) {
      return(NaN)
    }
    logS <- -c1 * q - (a / b) * (exp(b * q) - 1)
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

## quantile function
qgompzMakeNim <- nimbleFunction(
  run = function(p = double(0), a = double(0),
                 b = double(0), c1 = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    ## check inputs
    if(a < 0 | b < 0 | c1 < 0) {
      return(NaN)
    }
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    w0 <- (a / c1) * exp(a / c1) * ((1 - p)^(-b / c1))
    w0 <- RLW0(w0)
    q <- (a / (b * c1)) - (1 / c1) * log(1 - p) - (1 / b) *  w0
    return(q)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dgompzMakeNim = list(
    BUGSdist = "dgompzMakeNim(a, b, c1)",
    Rdist = "dgompzMakeNim(a, b, c1)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))
