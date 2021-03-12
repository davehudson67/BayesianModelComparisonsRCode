## create survival times
a <- 0.0000067
b <- 0.125
c1 <- 0.011

## simulate data
x <- seq(0, 120, length.out = 500)
sims <- rGompzMake(n, a, b, c1)
trued <- dGompzMake(x, a, b, c1)
trues <- pGompzMake(x, a, b, c1, lower.tail = FALSE)
trueh <- trued / trues
par(mfrow = c(1, 3))
hist(sims, freq = F, main = "PDF")
lines(x, trued, type = "l", col = "red")
plot(x, trueh, type = "l", main = "Hazard")
plot(x, trues, type = "l", main = "Survivor")
par(mfrow = c(1, 1))

## set up data
nind <- length(sims)
tD <- sims

while(length(which(is.infinite(tD))) > 0){
  inf_error <- which(is.infinite(tD))
  tD[inf_error] <- rGompzMake(length(inf_error), a, b, c1)
} 