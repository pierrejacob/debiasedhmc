library(debiasedhmc)
library(coda)
library(latex2exp)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()


### the target with constraints defined by linear inequalities
load("tmg.linear.samples.RData")
plot(samples, pch=".")
nsamples <- 1000
g <- qplot(samples[1:nsamples,1], samples[1:nsamples,2], geom = "point", alpha = 0.25)
g <- g + xlab(expression(x[1])) + ylab(expression(x[2])) + theme(legend.position = "none")
g
ggsave(filename = "tmg.dist1.png", plot = g, width = 9, height = 6)

load(file = "tmg.linear.cchains.RData")
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
hist(meetingtimes)
ghist1 <- qplot(x = meetingtimes, geom = "blank") + geom_histogram(aes(y = ..density..), binwidth = 1)
ghist1 <- ghist1 + xlab("meeting times")
ghist1
ggsave(filename = "tmg.meetingtimes1.pdf", plot = ghist1, width = 9, height = 6)

test_function <- function(x) x[1] > 3 & x[1] < 4
test_function <- function(x) x[2]

estim <- foreach(irep = 1:nrep, .combine = rbind) %dopar% {
  H_bar(cchains[[irep]], h = test_function, k = 10, K = K)
}
# confidence interval on the integral of interest
cat(colMeans(estim) - 2*sqrt(cov(estim) / nrep), colMeans(estim) + 2*sqrt(cov(estim) / nrep), "\n")
# MCMC estimate
mean(apply(samples[200:nrow(samples),], 1, test_function))


### the target with constraints defined by quadratic inequalities
load("tmg.quadratic.samples.RData")
plot(samples, pch=".")
nsamples <- 3000
g <- qplot(samples[1:nsamples,1], samples[1:nsamples,2], geom = "point", alpha = 0.25)
g <- g + xlab(expression(x[1])) + ylab(expression(x[2])) + theme(legend.position = "none")
g
ggsave(filename = "tmg.dist2.png", plot = g, width = 9, height = 6)


load(file = "tmg.quadratic.cchains.RData")
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
ghist2 <- qplot(x = meetingtimes, geom = "blank") + geom_histogram(aes(y = ..density..), binwidth = 1)
ghist2 <- ghist2 + xlab("meeting times")
ghist2
ggsave(filename = "tmg.meetingtimes2.pdf", plot = ghist1, width = 9, height = 6)

# test_function <- function(x) x[1] > 3 & x[1] < 4
# test_function <- function(x) x[2]
# estim <- foreach(irep = 1:nrep, .combine = rbind) %dopar% {
#   H_bar(cchains[[irep]], h = test_function, k = 10, K = K)
# }
# # confidence interval on the integral of interest
# cat(colMeans(estim) - 2*sqrt(cov(estim) / nrep), colMeans(estim) + 2*sqrt(cov(estim) / nrep), "\n")
# # MCMC estimate
# mean(apply(samples[200:nrow(samples),], 1, test_function))
#
