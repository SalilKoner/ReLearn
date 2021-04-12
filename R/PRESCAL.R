#' RESCAL with Projection using Alternate Least Squares
#'
#' Produce the the reduced rank factor \code{A} and the weight matrix \code{R} as a list for to each relation
#' @author Salil Koner, Dhrubajyoti Ghosh \cr Maintainer: Salil Koner
#' \email{skoner@@ncsu.edu}
#' @import Matrix
#' @import SparseM
#' @param X input list of length r , r being the number of relations.
#'          Each element of the list must be a n by n matrix. n being the number of entities
#' @param rank rank of the factorization
#' @param P Projection matrix
#' @param maxIter Number of iterations, default 500
#' @param conv Tolerance parameter for the convergence, default 1e-4
#' @param lam.A Tuning parameter for A
#' @param lam.R Tuning parameter for R
#' @param A0 Initial value for the n by rank matrix A
#' @param initA The initialization process of A : n by rank matrix. Can be either 'random' or 'EigenDecomp'. N
#'        Doesn't need to specified if A0 is provided. Default is 'EigenDecomp'
#' @param compute_fit Specifies whether the convergence is based on taking the difference of
#'                    the norm of X[[k]] and AR[[k]]A^T. If TRUE it computes the difference,
#'                     if n is very large, it is recommended to set FALSE.
#' @return The factorization A and a list of weights matrix R with length of the list as the length of X
#' @export
#' @examples
#' set.seed(12345)
#' n <- 4; p <- 2; d <- 2;
#' X1 <- sparseMatrix(i=c(1,3,4), j=c(3,2,4), x=1, dims=c(n,n))
#' X2 <- sparseMatrix(i=c(1,3,4,4), j=c(1,2,3,4), x=1, dims=c(n,n))
#' X  <- list(X1,X2)
#' result <- PRESCAL_ALS(X,p,P=matrix(rnorm(d*n),d,n), lam.A=1, lam.R=1, compute_fit=TRUE)
PRESCAL_ALS <- function(X, rank, P, ...){

  Args         <- list(...)
  ArgsProvided <- match.call()

  if (!(("X" %in% names(ArgsProvided)) & ("rank" %in% names(ArgsProvided)) & ("P" %in% names(ArgsProvided))) ){
    stop("Either X or the rank of the factorization or Projection matrix P is not provided")
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
  if (length(dim(P)) != 2){
    stop("P must be a matrix")
  }
  d              <- nrow(P)
  if (ncol(P) != n){
    stop("The number of columns of P must match the row/col of the frontal slices")
  }

  lam.A          <- ifelse("lam.A" %in% names(ArgsProvided), Args$lam.A, 1)
  lam.R          <- ifelse("lam.R" %in% names(ArgsProvided), Args$lam.R, 1)
  maxIter        <- ifelse("maxIter" %in% names(ArgsProvided), Args$maxIter, 500)
  initA          <- ifelse("initA" %in% names(ArgsProvided), Args$initA, "EigenDecomp")
  conv           <- ifelse("conv" %in% names(ArgsProvided), Args$conv, 1e-4)
  compute_fit    <- ifelse("compute_fit" %in% names(ArgsProvided), Args$compute_fit, FALSE)

  # print(lam.A); print(lam.R); print(maxIter) ; print(initA) ; print(conv) ; print(compute_fit);

  PX             <- lapply(1:r, function(k) as.matrix(P %*% X[[k]]) )
  if (compute_fit){
    normPX       <- 0
    for (k in 1:r){
      normPX     <- normPX + sum(PX[[k]]^2)
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

  normP <- (max(svd(P)$d)^2)

  fit <- fitchange <- fitold <- 0
  for(itr in 1:maxIter){

    if (itr %% 10==0){
      print(paste0("iteration:", itr))
    }


    fitold     <- fit

    AtA        <- crossprod(A)
    PA         <- P %*% A

    # Update R
    svdA      <- svd(A) ; svdPA <- svd(PA);
    evs       <- kronecker(svdA$d, svdPA$d);
    D         <- matrix(evs/(evs^2 + lam.R), min(d,p), p)
    Rs         <- list()
    for (k in 1:r){
      R        <- D * (crossprod(svdPA$u, PX[[k]]) %*% svdA$u)
      Rs[[k]]  <- (svdPA$v)%*%tcrossprod(R, svdA$v)
    }


    left       <- matrix(0, p, p)
    right      <- matrix(0, p, p)
    for (k in 1:r){
      left     <- left  + (Rs[[k]] %*% AtA %*% t(Rs[[k]]))
      right    <- right + crossprod(PA %*% Rs[[k]])
    }

    norm.left  <- max(svd(left)$d)*normP ; norm.right   <- max(svd(right)$d) ;
    max.step   <- 1/(norm.left + norm.right) ; step.size  <- 0.5*max.step
    # norm       <- max(svd(left + right)$d) ; # print(norm)
    # max.step   <- 1/norm ; step.size    <- max.step/2

    grad_Plus  <- t(P) %*% (PA %*% left) + A %*% right + lam.A*A
    left       <- matrix(0, d, p)
    right      <- matrix(0, n, p)
    for (k in 1:r){
      left     <- left  + PX[[k]] %*% tcrossprod(A, Rs[[k]])
      right    <- right + crossprod(PX[[k]], PA %*% Rs[[k]])
    }
    grad_Minus <- t(P) %*% left + right
    grad       <- grad_Plus - grad_Minus
    # print(grad)
    # print(max(svd(grad)$d))
    A          <- A - step.size*grad

    if (compute_fit){
      normdiff    <- 0
      for (k in 1:r){
        PARAt     <- PA %*% Rs[[k]] %*% t(A)
        normdiff  <- normdiff + (norm(PX[[k]] - PARAt, type="F"))^2
      }
      fit      <- 1 - normdiff/normPX
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
