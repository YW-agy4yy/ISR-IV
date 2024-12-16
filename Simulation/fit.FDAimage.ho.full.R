fit.FDAimage.ho.full <- function (Y, X, Z, V, Tr, d = 5, r = 1, lambda = 0.1) 
{
  if (!is.matrix(Y)) {
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  if (!is.matrix(X)) {
    warning("The explanatory variable, X, should be a matrix.")
    X <- as.matrix(X)
  }
  np <- ncol(X)
  if (!is.matrix(Z)) {
    warning("The coordinates of each pixel/voxel, Z, should be a matrix.")
    Z <- as.matrix(Z)
  }
  lambda <- as.matrix(lambda)
  if (nrow(lambda) > 1) {
    warning("The tuning parameter, lambda, should be a scalar. Instead, the default 0.1 is used.")
    lambda <- as.matrix(0.1)
  }
  this.call <- match.call()
  Ball <- basis(V, Tr, d, r, Z)
  K <- Ball$K
  Q2 <- Ball$Q2
  B <- Ball$B
  ind.inside <- Ball$Ind.inside
  npix <- length(ind.inside)
  tria.all <- Ball$tria.all
  n = nrow(Y)
  nx = ncol(X)
  npix = ncol(Y)
  BQ2 = B %*% Q2
  BQ2 = as.matrix(BQ2)
  K = as.matrix(K)
  J = ncol(BQ2)
  WW = kronecker(crossprod(X), crossprod(BQ2))
  rhs = rowSums(sapply(1:n, function(iter) as.matrix(kronecker(X[iter, 
  ], crossprod(BQ2, Y[iter, ])))))
  P = as.matrix(crossprod(Q2, K) %*% Q2)
  lambda = as.matrix(lambda)
  nlam = nrow(lambda)
  gamma_all = list(nlam)
  beta_all = list(nlam)
  sse_all = c()
  for (il in 1:nlam) {
    Lam = diag(rep(lambda[il], nx))
    Dlam = as.matrix(kronecker(Lam, P))
    lhs = WW + Dlam
    theta = solve(lhs, rhs)
    theta.mtx = matrix(theta, J, nx)
    gamma = Q2 %*% theta.mtx
    gamma_all[[il]] = gamma
    beta = BQ2 %*% theta.mtx
    beta_all[[il]] = beta
    Yhat = tcrossprod(X, beta)
    ssei = apply((Y - Yhat)^2, 1, sum)
    sse = sum(ssei)
    sse_all = c(sse_all, sse)
  }
  list(beta = beta_all, gamma = gamma_all, sse = sse, theta.mtx = theta.mtx)
}