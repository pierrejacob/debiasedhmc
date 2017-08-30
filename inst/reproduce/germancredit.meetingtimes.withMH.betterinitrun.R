library(debiasedhmc)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()
registerDoParallel(cores = detectCores())
##
get_hmc_kernel <- function(logtarget, gradlogtarget, stepsize, nsteps, dimension){
  # leap frog integrator
  # note that some papers use the notation U for - logtarget, so that there are minus signs everywhere
  leapfrog <- function(x, v){
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }
  # One step of HMC
  kernel <- function(chain_state, iteration){
    current_v <- rnorm(dimension) # velocity or momentum
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x
    # if we were a bit smarter, we would save logtarget(chain_state) so as
    # to not re-evaluate it at every step
    accept_ratio <- logtarget(proposed_x) - logtarget(chain_state)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
    accept <- FALSE
    if (log(runif(1)) < accept_ratio){
      chain_state <- proposed_x
      current_v <- proposed_v
      accept <- TRUE
    } else {
    }
    return(list(chain_state = chain_state, accept = accept))
  }
  # One step of HMC
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    current_v <- rnorm(dimension) # velocity or momentum, shared by both chains
    leapfrog_result1 <- leapfrog(chain_state1, current_v)
    leapfrog_result2 <- leapfrog(chain_state2, current_v)
    proposed_v1 <- - leapfrog_result1$v
    proposed_x1 <- leapfrog_result1$x
    proposed_v2 <- - leapfrog_result2$v
    proposed_x2 <- leapfrog_result2$x
    # if we were a bit smarter, we would save logtarget(current_x) so as
    # to not re-evaluate it at every step
    accept_ratio1 <- logtarget(proposed_x1) - logtarget(chain_state1)
    accept_ratio2 <- logtarget(proposed_x2) - logtarget(chain_state2)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio1 <- accept_ratio1 + sum(current_v^2) / 2 - sum(proposed_v1^2) / 2
    accept_ratio2 <- accept_ratio2 + sum(current_v^2) / 2 - sum(proposed_v2^2) / 2
    logu <- log(runif(1)) # shared by both chains
    if (logu < accept_ratio1){
      chain_state1 <- proposed_x1
      current_v1 <- proposed_v1
    } else {
    }
    if (logu < accept_ratio2){
      chain_state2 <- proposed_x2
      current_v2 <- proposed_v2
    } else {
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}

germancredit <- read.table(system.file("", "germancredit.txt", package = "debiasedhmc"))

X <- germancredit[,1:24]
X <- scale(X)
Y <- germancredit[,25]
Y <- Y - 1
n <- nrow(X)
p <- ncol(X)

Xinter <- matrix(nrow = n, ncol = p*(p-1) / 2)
index <- 1
for (j in 1:(p-1)){
  for (jprime in (j+1):p){
    Xinter[,index] <- X[,j] * X[,jprime]
    index <- index + 1
  }
}
Xinter <- scale(Xinter)
X <- cbind(X, Xinter)
p <- ncol(X)
lambda <- 0.01
dimension <- p + 2
# create initial distribution from approx of target
# function to compute distance between coupled chains
distance_ <- function(cchain){
  nsteps <- nrow(cchain$samples1) - 1
  return(sapply(1:nsteps, function(index) sqrt(sum((cchain$samples1[1+index,] - cchain$samples2[index,])^2))))
}
##
load(file = "germancredit.withMH.cchains.RData")
distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  ds <- distance_(cchains[[irep]])
  if (length(ds) < max_iterations){
    ds <- c(ds, rep(0, max_iterations - length(ds)))
  }
  ds
}

matplot(log(t(distances_)), type = "l", xlim = c(1,500), col = "black")
k <- 100
mean_estimators <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  H_bar(cchains[[irep]], k = k, K = K)
}
square_estimators <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  H_bar(cchains[[irep]], h = function(x) x^2, k = k, K = K)
}
post_mean <- colMeans(mean_estimators)
post_var <- colMeans(square_estimators) - post_mean^2
# need to make sure that post_vars only have positive entries
post_var <- diag(post_var, dimension, dimension)
summary(diag(post_var))

# all_samples <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
#   cchains[[irep]]$samples1[(K-100+2):(K+1),]
# }
# post_mean <- colMeans(all_samples)
# post_var <- cov(all_samples)
# post_var <- diag(diag(post_var))
rinit <- function() fast_rmvnorm(1, post_mean, post_var)[1,]
##
logtarget <- function(x) debiasedhmc:::logistic_logtarget_c(x, Y, X, lambda)
gradlogtarget <- function(x) debiasedhmc:::logistic_gradlogtarget_c(x, Y, X, lambda)

nsteps <- 20
effectiveTime <- 0.1
cat("effective time:", effectiveTime, "\n")
stepsize <- effectiveTime / nsteps
hmc_kernel <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

omega <- 1/20
Sigma_proposal <- 1e-10 * diag(1, dimension, dimension)
mh_kernel <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# mixture kernels
mixture_kernel <- function(chain_state, iter){
  if (runif(1) < omega){
    return(mh_kernel$kernel(chain_state, iter))
  } else {
    return(hmc_kernel$kernel(chain_state, iter))
  }
}
#
mixture_coupled_kernel <- function(chain_state1, chain_state2, iter){
  if (runif(1) < omega){
    return(mh_kernel$coupled_kernel(chain_state1, chain_state2, iter))
  } else {
    return(hmc_kernel$coupled_kernel(chain_state1, chain_state2, iter))
  }
}


nrep <- 100
max_iterations <- 1000
K <- 1000
cchains <- foreach(irep = 1:nrep) %dorng% {
  coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K, max_iterations = max_iterations)
}
save(K, max_iterations, nrep, effectiveTime, nsteps, stepsize, cchains, file = "germancredit.withMH.cchains.betterinit.RData")
load(file = "germancredit.withMH.cchains.betterinit.RData")

nrep <- 1000
K <- 1000
cchains <- foreach(irep = 1:nrep) %dorng% {
  coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K)
}
save(K, nrep, effectiveTime, nsteps, stepsize, cchains, file = "germancredit.withMH.cchains.betterinit.1000.RData")
# load(file = "germancredit.withMH.cchains.betterinit.1000.RData")

### The following function returns H_ell for ell = k ... K, in order O(K) operations
#
# H_all <- function(c_chains, h = function(x) x, k = 0, K = 1){
#   maxiter <- c_chains$iteration
#   if (k > maxiter){
#     print("error: k has to be less than the horizon of the coupled chains")
#     return(NULL)
#   }
#   if (K > maxiter){
#     print("error: K has to be less than the horizon of the coupled chains")
#     return(NULL)
#   }
#   # test the dimension of h(X)
#   h_of_chain <- apply(X = c_chains$samples1[(k+1):(K+1),,drop=F], MARGIN = 1, FUN = h)
#   if (is.null(dim(h_of_chain))){
#     h_of_chain <- matrix(h_of_chain, ncol = 1)
#   } else {
#     h_of_chain <- t(h_of_chain)
#   }
#   H_all <- matrix(nrow = nrow(h_of_chain), ncol = ncol(h_of_chain))
#   sumDeltas <- rep(0, ncol(h_of_chain))
#   H_all[nrow(H_all),] <- h(c_chains$samples1[K + 1,])
#   for (t in rev(k:(K-1))){
#     sumDeltas <- sumDeltas + h(c_chains$samples1[t + 2,]) - h(c_chains$samples2[t+1,])
#     H_all[t-k+1,] <- sumDeltas + h(c_chains$samples1[t + 1,])
#   }
#   return(H_all)
# }
