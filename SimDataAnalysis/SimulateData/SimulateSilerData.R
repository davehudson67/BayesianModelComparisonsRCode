## create survival times
a1 <- exp(-2)
a2 <- exp(-11)
b1 <- 1.2
b2 <- 0.11
c <- exp(-5)

## simulate data
x <- seq(0, 120, length.out = 500)
sims <- rSiler(n, a1, a2, b1, b2, c)
trued <- dSiler(x, a1, a2, b1, b2, c)
trues <- pSiler(x, a1, a2, b1, b2, c, lower.tail = FALSE)
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