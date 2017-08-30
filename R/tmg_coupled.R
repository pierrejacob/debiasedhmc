### This file is taken/modified from the tmg package on CRAN, by Ari Pakman
#'@export
rtmg_coupled_R <- function(n, M, r, initial1, initial2, f = NULL, g = NULL, q = NULL){
  burn.in <- 0

  d = nrow(M)      # dimension of target space
  if (ncol(M)!=d){
    cat("Error: M must be a square matrix.")
    return()
  }

  if (length(initial1)!=d){
    cat("Error: wrong length for initial value vector.")
    return()
  }

  # symmetrize M and verify that it is positive definite
  M = (M + t(M))/2;
  eigs=eigen(M, symmetric=T, only.values=T)$values
  if (any(eigs<=0)){
    cat("Error: M must be positive definite.")
    return()
  }


  # we change variables to the canonical frame, sample by calling the C++ code and transform back to the original frame

  R = chol(M)
  Mir = solve(M,r)
  Ri = solve(R)
  initial12= as.vector(R%*%initial1) - as.vector(r%*%Ri)
  initial22= as.vector(R%*%initial2) - as.vector(r%*%Ri)

  if (is.null(f) & is.null(g)){
    numlin=0
    f2=NULL
    g2=NULL
  } else   # if there are linear constraints
  { if (is.matrix(f) & is.vector(g)){
    # verify linear constraints sizes
    numlin = nrow(f);
    if (length(g) != numlin | ncol(f) != d ){
      cat("Error: inconsistent linear constraints. f must be an m-by-d matrix and g an m-dimensional vector.")
      return()
    }

    # verify initial value satisfies linear constraints
    if (any(f%*%initial1 + g <= 0)) {
      cat("Error: initial point violates linear constraints.")
      return()
    }

    # map linear constraints to the canonical frame
    f2 = f%*%Ri
    g2 = as.vector(f%*%Mir + g)

  } # if (is.matrix(f) & is.vector(g))
    else {
      cat("Error: for linear constraints, f must be a matrix and g a vector.\n")
      return()
    }}


  if (is.null(q)){
    numquad=0
    quads = NULL
  } else    # if there are quadratic constraints
  { if (is.list(q)){
    # verify that the elements of the quadratic constraints are lists of length  3
    ll=lapply(q,length)
    nc =c(do.call("cbind",ll))
    if (any(nc!=3)) {
      cat("Error: each element in q must be a list of length 3.\n");
      return()
    }

    numquad = length(q);
    quads = matrix( nrow=numquad*(d+2), ncol=d)

    for (i in 1:numquad){
      qci = q[[i]];
      t1 = initial1 %*% qci[[1]] %*% initial1 + qci[[2]] %*% initial1 + qci[[3]];
      t2 = initial2 %*% qci[[1]] %*% initial2 + qci[[2]] %*% initial2 + qci[[3]];
      if (t1 <=0 || t2 <= 0){
        cat("Error: initial point violates quadratic constraints. \n")
        return()
      } else {

        # map quadratic constraints to the canonical frame
        A= qci[[1]]
        B= qci[[2]]
        C= qci[[3]]
        quads[ ((i-1)*(d+2)+1): ((i-1)*(d+2)+d)  , ] = t(Ri)%*% A %*% Ri
        quads[ ((i-1)*(d+2)+d+1) ,   ] = 2*t(Mir) %*% A %*% Ri + t(B) %*% Ri
        C = C + t(Mir)%*%A%*% Mir + t(B)%*% Mir
        quads[ i*(d+2) ,  ] = c(C, rep(0,d-1))
      }
    } #for (i in 1:numquad)
  } #if (is.list(q))
    else {
      cat("Error: for quadratic constraints, q must be a list.\n")
      return()
    }}

  seed = sample(1:10000000,1)
  if (is.null(f2)) f2 <- matrix()
  if (is.null(g2)) g2 <- vector()
  if (is.null(quads)) quads <- matrix()

  samples = rtmg_coupled(n+burn.in, seed, initial12, initial22, numlin, f2, g2, numquad, quads)
  samples$samples1 = samples$samples1%*%t(Ri) + matrix(rep(Mir,n+1),nrow=n+1, byrow=T)
  samples$samples2 = samples$samples2%*%t(Ri) + matrix(rep(Mir,n),nrow=n, byrow=T)
  return(samples)
}
