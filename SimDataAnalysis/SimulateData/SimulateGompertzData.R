## create survival times
a <- 0.121
b <- 0.105

## simulate data
x <- seq(0, 120, length.out = 500)
sims <- rGompertz(n , a, b)
trued <- dGompertz(x, a, b)
trues <- pGompertz(x, a, b, lower.tail = FALSE)
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