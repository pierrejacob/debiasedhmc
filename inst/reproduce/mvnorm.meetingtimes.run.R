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
effectiveTime <- pi/2
nsteps <- 20
stepsize <- effectiveTime / nsteps
hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)

nrep <- 100
max_iterations <- 10000
cchains <- foreach(irep = 1:nrep) %dorng% {
  coupled_chains(hmc_kernel$kernel, hmc_kernel$coupled_kernel, rinit, max_iterations = max_iterations)
}
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
save(max_iterations, nrep, effectiveTime, nsteps, stepsize, meetingtimes, cchains, file = "mvnorm.cchains.RData")


# load(file = "mvnorm.cchains.RData")
# summary(meetingtimes)
#
# distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
#   ds <- distance_(cchains[[irep]])
#   if (length(ds) < max_iterations){
#     ds <- c(ds, rep(0, max_iterations - length(ds)))
#   }
#   ds
# }
# matplot(log(t(distances_)), type = "l", xlim = c(1,300), col = "black")
