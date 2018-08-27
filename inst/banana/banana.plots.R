rm(list = ls())
library(debiasedhmc)
library(doParallel)
library(doRNG)
library(ggplot2)
library(gridExtra)
setmytheme()
registerDoParallel(cores = detectCores()-2)
set.seed(1)

load("banana.synchr.meetings.RData")

meetings.synchronous <- sapply(results.df, function(x) x$meetingtime)

load("banana.contr.meetings.RData")
meetings.contractive <- sapply(results.df, function(x) x$meetingtime)

#
summary(meetings.synchronous)
summary(meetings.contractive)

gsync <- qplot(x = meetings.synchronous, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("Meeting times") + ylab("Density")
ggsave(filename = "banana.synchr.meetings.pdf", plot = gsync)

gcontr <- qplot(x = meetings.contractive, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("Meeting times") + ylab("Density")
ggsave(filename = "banana.contr.meetings.pdf", plot = gcontr)


