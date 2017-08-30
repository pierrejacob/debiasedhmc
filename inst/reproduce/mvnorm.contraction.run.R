library(debiasedhmc)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
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
#
effectiveTimes <- c(pi/4, 3*pi/8, 7*pi/16, 2*pi/4, 9*pi/16, 5*pi/8, 3*pi/4, pi, 5*pi/4, 6*pi/4)

neff <- length(effectiveTimes)
nsteps <- 20
nrep <- 5
maxiterations <- 100
#
distance.df <- data.frame()
for (i in 1:neff){
  effectiveTime <- effectiveTimes[i]
  cat("effective time:", effectiveTime, "\n")
  stepsize <- effectiveTime / nsteps
  hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)
  cchains <- foreach(irep = 1:nrep) %dorng% {
    coupled_chains(hmc_kernel$kernel, hmc_kernel$coupled_kernel, rinit, max_iterations = maxiterations)
  }
  distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
    ds <- distance_(cchains[[irep]])
    if (length(ds) < maxiterations){
      ds <- c(ds, rep(0, maxiterations - length(ds)))
    }
    ds
  }
  distance.df <- rbind(distance.df, data.frame(irep = 1:nrep, effectiveTime = rep(effectiveTime, nrep),
             stepsize = rep(stepsize, nrep),
             nsteps = rep(nsteps, nrep),
             distance = distances_[,maxiterations]))
  save(maxiterations, nrep, nsteps, effectiveTimes, distance.df, file = "mvnorm.contraction.RData")
}

# distance.df
# load(file = "mvnorm.contraction.RData")
# g <- ggplot(distance.df %>% filter(distance > 1e-30), aes(x = effectiveTime, y = distance)) + geom_point()
# g <- g + scale_y_log10(breaks = c(1e-20, 1e-10, 1), limits = c(1e-30, 10))
# g

