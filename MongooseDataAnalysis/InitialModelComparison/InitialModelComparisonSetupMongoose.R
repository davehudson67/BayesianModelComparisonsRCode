## load libraries
library(tidyverse)
library(survminer)
library(survival)
rm(list=ls())

## load data
mong <- readRDS("../Data/mongoose.rds")

## set up data
tU <- mong$death_day
tB <- mong$birth_day
tL <- mong$last_seen

## check last seen against death
table(tL >= tU)
summary(mong$birth_day[tL >= tU])
summary(mong$death_day[tL < tU])

## set death time for censored individuals
#tU[tU == 0] <- NA

## set up censoring vector (1 = interval, 2 = right)
censored <- mong$censored

## summaries
stopifnot(all(tL[!is.na(tU)] <= tU[!is.na(tU)]))
stopifnot(all(tB <= tL)) #some individuals died on the same day they were born

## normalise to survival times
tU <- tU - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tU)
colnames(cint) <- NULL
cint[censored == 2, 2] <- cint[censored == 2, 1] 
cint[censored == 2, 1] <- 0

## check censoring times
summary(apply(cbind(tL, tU), 1, diff)[censored == 1])
summary(apply(cbind(tL, tU), 1, diff)[censored == 2])

## set up latent death times
tD <- rep(NA, length(tU))

## set up nind
nind <- nrow(cint)

## set up zL/zU input matrix
zL <- tL
zU <- tU

## set up zL/zU input matrix
zL <- tL
zU <- tU

## produce K-M plot (assuming that interval-censored individuals
## die at the end time of the interval)
mong_km_dat <- select(mong, tD = death_day, tL = last_seen, tB = birth_day, cens = censored) %>%
  mutate(tD = ifelse(cens == 2, tL, tD)) %>%
  mutate(tD = tD - tB) %>%
  select(-tL, -tB) %>%
  mutate(cens = 1 - (cens - 1))

## plot non-sex specific KM
temp_fit <- surv_fit(Surv(tD, cens) ~ 1, data = mong_km_dat)
km_plot <- ggsurvplot(temp_fit, data = mong_km_dat, conf.int = TRUE)
km_plot

## set and sample some random seeds for different models
set.seed(42)
seeds <- round(runif(5, 0, 100000000))

## create output folder
dir.create("outputs")

## save data
rm(temp_fit, tL, tU, mong, tB)
save.image("mong.RData")