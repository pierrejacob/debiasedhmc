rm(list = ls())
library(ggplot2)
library(tictoc)
library(debiasedhmc)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
setmytheme()
set.seed(1)


## Plot for RWMH
load(file = "scaling.rwmh.meetings.v3.RData")

ggplot(results.df, aes(x = dimension, y = mean_meetingtime, group = constant, linetype = factor(constant))) + geom_line() +
  scale_linetype(name = "C =") + ylab("Average meeting time") + xlab("Dimension") + ylim(0,3000)


## Plot for MALA
load(file = "scaling.mala.meetings.v3.RData")

ggplot(results.df, aes(x = dimension, y = mean_meetingtime, group = constant, linetype = factor(constant))) + geom_line() +
  scale_linetype(name = "C =") + ylab("Average meeting time") + xlab("Dimension")

## Plot for HMC
load(file = "scaling.hmc.meetings.v3.RData")

ggplot(results.df, aes(x = dimension, y = mean_meetingtime, group = constant, linetype = factor(constant))) + geom_line() +
  scale_linetype(name = "C =") + ylab("Average meeting time") + xlab("Dimension")
