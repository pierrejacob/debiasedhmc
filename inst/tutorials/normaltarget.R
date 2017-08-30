library(debiasedhmc)
rm(list = ls())
set.seed(18)
setmytheme()
registerDoParallel(cores = detectCores())
# define target of dimension p
dimension <- 250
mean_target <- rep(0, dimension)
Sigma_target <- diag(1, dimension, dimension)
for (i in 1:dimension){
  for (j in 1:dimension){
    Sigma_target[i,j] <- exp(-abs(i-j))
  }
}
target <- get_mvnormal(dimension, mean_target, Sigma_target)
# initial distribution of chains
rinit <- function() fast_rmvnorm(1, mean_target, Sigma_target)
# function to compute distance between coupled chains
distance_ <- function(cchain){
  nsteps <- nrow(cchain$samples1) - 1
  return(sapply(1:nsteps, function(index) sqrt(sum((cchain$samples1[1+index,] - cchain$samples2[index,])^2))))
}

# first, let's generate pairs of long Hamiltonian trajectories
# with common initial velocity

# leap-frog step size
stepsize <- 0.05
trajectorylength <- 10
nsteps <- floor(trajectorylength/stepsize)
# load HMC kenrel
hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)
# initial positions
x1 <- rinit()
x2 <- rinit()

#
result <- hmc_kernel$coupled_kernel(x1, x2)
# trajectory of particle 1
xtraj1 <- result$xtraj1
# trajectory of particle 2
xtraj2 <- result$xtraj2
xtraj1.df <- data.frame(xtraj1)
xtraj2.df <- data.frame(xtraj2)
g <- ggplot(xtraj1.df, aes(x = X1, y = X2, colour = "trajectory 1")) + geom_path() + theme(legend.position = "none")
g <- g + geom_path(data = xtraj2.df, aes(colour = "trajectory 2"))
g <- g + xlab(expression(x[1])) + ylab(expression(x[2]))
print(g)
# nothing in particular happens!
# we can compute the distance between the two trajectories
distances <- sapply(1:(nsteps+1), FUN = function(index) sqrt(sum((xtraj1[index,] - xtraj2[index,])^2)))
gdistance <- qplot(x = 0:nsteps, y = distances, geom = "line") + xlab("iteration") + ylab("distance")
print(gdistance)
# it fluctuates, but note that it starts by decreasing

### Now we look at what happens when we iterate short trajectories,
### each time redrawing the velocity...

x1 <- rinit()
x2 <- rinit()
trajectorylength <- 1
nsteps <- floor(trajectorylength/stepsize)
hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)

ntrajectories <- 20
xtraj1 <- matrix(nrow = ntrajectories * (nsteps+1), ncol = dimension)
xtraj2 <- matrix(nrow = ntrajectories * (nsteps+1), ncol = dimension)
for (itraj in 1:ntrajectories){
  result <- hmc_kernel$coupled_kernel(x1, x2)
  x1 <- result$chain_state1
  x2 <- result$chain_state2
  xtraj1[((itraj-1)*(nsteps+1)+1):(itraj*(nsteps+1)),] <- result$xtraj1
  xtraj2[((itraj-1)*(nsteps+1)+1):(itraj*(nsteps+1)),] <- result$xtraj2
}

xtraj1.df <- data.frame(xtraj1)
xtraj2.df <- data.frame(xtraj2)
xtraj1.df$iteration <- rep(1:ntrajectories, each = (nsteps+1))
xtraj1.df$step <- rep(0:nsteps, times = ntrajectories)
xtraj2.df$iteration <- rep(1:ntrajectories, each = (nsteps+1))
xtraj2.df$step <- rep(0:nsteps, times = ntrajectories)

g <- ggplot(xtraj1.df, aes(x = X1, y = X2, colour = "trajectory 1")) + geom_path() + theme(legend.position = "none")
g <- g + geom_path(data = xtraj2.df, aes(colour = "trajectory 2"))
g <- g + xlab(expression(x[1])) + ylab(expression(x[2]))
g <- g + geom_point(aes(x = xtraj1[1,1], y = xtraj1[1,2], colour = "trajectory 1"), size = 6)
g <- g + geom_point(aes(x = xtraj2[1,1], y = xtraj2[1,2], colour = "trajectory 2"), size = 6)
print(g)
# the two trajectories get closer to one another!

distances <- sapply(1:nrow(xtraj1), FUN = function(index) sqrt(sum((xtraj1[index,] - xtraj2[index,])^2)))
gdistance <- qplot(x = 1:nrow(xtraj1) / nsteps, y = distances, geom = "line") + xlab("iteration") + ylab("distance")
print(gdistance)

