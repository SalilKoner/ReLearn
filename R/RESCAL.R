#' RESCAL algorithm using Alternating Least Squares
#'
#' Produce the the reduced rank factor \code{A} and the weight matrix \code{R} as a list for to each relation
#' @author Salil Koner, Dhrubajyoti Ghosh \cr Maintainer: Salil Koner
#' \email{skoner@@ncsu.edu}
#' @import Matrix
#' @import SparseM
#' @param X input list of length r , r being the number of relations.
#' Each element of the list must be a n by n matrix. n being the number of entities
#' @param rank rank of the factorization
#' @param P Projection matrix
#' @param maxIter Number of iterations, default 500
#' @param conv Tolerance parameter for the convergence, default 1e-4
#' @param lam.A Tuning parameter for A
#' @param lam.R Tuning parameter for R
#' @param A0 Initial value for the n by rank matrix A
#' @param initA The initialization process of A : n by rank matrix. Can be either 'random' or 'EigenDecomp'. N
#' Doesn't need to specified if A0 is provided. Default is 'EigenDecomp'
#' @param compute_fit Specifies whether the convergence is based on taking the difference of
#' the norm of X[[k]] and AR[[k]]A^T. If TRUE it computes the difference,
#' if n is very large, it is recommended to set FALSE.
#' @return The factorization A and a list of weights matrix R with length of the list as the length of X
#' @export
#' @examples
#' set.seed(12345)
#' n <- 4; p <- 2; d <- 2;
#' X1 <- sparseMatrix(i=c(1,3,4), j=c(3,2,4), x=1, dims=c(n,n))
#' X2 <- sparseMatrix(i=c(1,3,4,4), j=c(1,2,3,4), x=1, dims=c(n,n))
#' X  <- list(X1,X2)
#' result <- RESCAL_ALS(X,p,lam.A=1, lam.R=1, compute_fit=TRUE)
#' @references Nickel, M., Tresp, V., & Kriegel, H. P. (2011, January)
#' \emph{A three-way model for collective learning on multi-relational data. In Icml} \cr
#' \url{http://www.icml-2011.org/papers/438_icmlpaper.pdf}.
RESCAL_ALS <- function(X, rank, ...){

  Args         <- list(...)
  # print(Args)

  # print(compute_fit)
  ArgsProvided <- match.call()


  if (!(("X" %in% names(ArgsProvided)) & ("rank" %in% names(ArgsProvided))) ){
    stop("Either X or the rank of the factorization is not provided")
  }

  r        <- length(X)

  for (k in 1:r){
    if (length(dim(X[[k]])) != 2){
      stop(paste0(toOrdinal::toOrdinal(k), " frontal slice is not a matrix"))
    }

    if (nrow(X[[k]]) != ncol(X[[k]])){
      stop(paste0(toOrdinal::toOrdinal(k), " frontal slice is not a square matrix"))
    }

    if (any(dim(X[[k]]) != dim(X[[1]]))){
      stop("Frontal slices of X must be all of same shape")
    }

    X[[k]] <- as(X[[k]], "dgCMatrix")
  }

  n        <- nrow(X[[1]]) ; p <- rank ;

  if (rank > n){
    stop("The rank of factorization can not be more than row/col dimension of slices")
  }

  lam.A          <- ifelse("lam.A" %in% names(ArgsProvided), Args$lam.A, 1)
  lam.R          <- ifelse("lam.R" %in% names(ArgsProvided), Args$lam.R, 1)
  maxIter        <- ifelse("maxIter" %in% names(ArgsProvided), Args$maxIter, 500)
  initA          <- ifelse("initA" %in% names(ArgsProvided), Args$initA, "EigenDecomp")
  conv           <- ifelse("conv" %in% names(ArgsProvided), Args$conv, 1e-4)
  compute_fit    <- ifelse("compute_fit" %in% names(ArgsProvided), Args$compute_fit, TRUE)



  if (compute_fit){
    normX        <- 0
    for (k in 1:r){
      normX      <- normX + sum(X[[k]]^2)
    }
  }

  if ("A0" %in% names(ArgsProvided)){
    A            <- A0
  }else{
    if (initA == "EigenDecomp"){
      X.add      <- Reduce("+", X) + Reduce("+", lapply(X, t))
      eigX.add   <- eigen(X.add)
      A          <- eigX.add$vectors
      A          <- A[,1:p]
    } else if (initA == "random"){
      A          <- matrix(rnorm(n*p), n, p)
    } else{
      stop("Method to initialize A can be only EigenDecomp or random, OR
            provide a initial value of A as an argument A0")
    }
  }

  fit <- fitchange <- fitold <- 0
  for(itr in 1:maxIter){

    if (itr %% 10==0){
      print(paste0("iteration:", itr))
    }

    fitold     <- fit

    AtA        <- crossprod(A)

    # Update R
    svdA       <- svd(A)
    evs        <- kronecker(svdA$d, svdA$d)
    D          <- matrix(evs/(evs^2 + lam.R), p, p)

    Rs         <- list()
    for (k in 1:r){
      R        <- D * (crossprod(svdA$u, X[[k]]) %*% svdA$u)
      Rs[[k]]  <- (svdA$v) %*% tcrossprod(R, svdA$v)
    }


    left       <- lam.A*diag(p)
    right      <- matrix(0, n, p)
    for (k in 1:r){
      left     <- left + Rs[[k]] %*% AtA %*% t(Rs[[k]])
      left     <- left + t(Rs[[k]]) %*% AtA %*% Rs[[k]]
      right    <- right + X[[k]] %*% tcrossprod(A, Rs[[k]])
      right    <- right + crossprod(X[[k]] , A %*% Rs[[k]])
    }


    A          <- t(solve(left, t(right)))

    if (compute_fit){
      normdiff    <- 0
      for (k in 1:r){
        ARAt      <- A %*% Rs[[k]] %*% t(A)
        normdiff  <- normdiff + (norm(X[[k]] - ARAt, type="F"))^2
      }
      fit      <- 1 - normdiff/normX
    }
    else{
      fit      <- itr
    }

    fitchange  <- abs(fit-fitold)
    if (fitchange < conv){
      break
    }

  }

  list("A"=A, "R"=Rs, "fit"=fit, "fitchange"=fitchange, "iter"=itr)

}

