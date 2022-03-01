#' Negative Laplacian of 2-D Radial Basis Function
#' @param kerngen kernel generator
#'
#' @export
neglap_rbf2D <- function(X, Y, kerngen){
  ep <- kerngen$ep
  drbf <- kerngen$drbf
  d2rbf <- kerngen$d2rbf
  r <- calc_distXYC(X,Y)
  A <- -d2rbf(ep,r)*(r^2) - 2*drbf(ep,r)
}
#' Build Gram matrix for divergence-free normal projection
#' @export
build_DFprojGram <- function(Xe, X, Xb, Nb, kerngen){
  A <- neglap_rbf2D(Xe, X, kerngen)
  B <- DF_rbf2D(Xe, Xb, kerngen)

  Dn <- rbind(Matrix::Diagonal(x = Nb[,1]),
              Matrix::Diagonal(x = Nb[,2]))
  B <- B %*% Dn

  C <- DF_rbf2D(Xb, Xb, kerngen);
  C <- Matrix::t(Dn) %*% C %*% Dn;
  return(list(A=A, B=B, C=C))
}

#' @export
DF_rbf2D <- function(X, Y, kerngen){
  ep <- kerngen$ep
  drbf <- kerngen$drbf
  d2rbf <- kerngen$d2rbf

  N <- nrow(X)
  M <- nrow(Y)

  r <- calc_distXYC(X, Y)

  D1 <- drbf(ep, r)
  D2 <- d2rbf(ep,r)


  rm(list = "r")
  xdiff <- outer(X[,1], Y[,1], '-')
  ydiff <- outer(X[,2], Y[,2], '-')

  A <- matrix(0, nrow =2*N, ncol =2*M)
  A[1:N, 1:M] <- -D2*ydiff^2 - D1
  A[1:N, (M+1):(2*M)] <- D2*xdiff*ydiff
  rm(list="ydiff")

  A[(N+1):(2*N), 1:M] <- A[1:N, (M+1):(2*M)]
  A[(N+1):(2*N),(M+1):(2*M)] <- -D2*xdiff^2 - D1
  return(A)
}

#' @export
CF_rbf2D <- function(X, Y, kerngen){
  ep <- kerngen$ep
  drbf <- kerngen$drbf
  d2rbf <- kerngen$d2rbf

  N <- nrow(X)
  M <- nrow(Y)

  r <- calc_distXYC(X, Y)

  D1 <- drbf(ep,r)
  D2 <- d2rbf(ep, r)

  rm(list = "r")
  xdiff <- outer(X[,1], Y[,1], '-')
  ydiff <- outer(X[,2], Y[,2], '-')

  A <- matrix(0, nrow =2*N, ncol =2*M)
  A[1:N, 1:M] <- -D2*xdiff^2 - D1
  A[1:N, (M+1):(2*M)] <- -D2*xdiff*ydiff
  rm(list="xdiff")

  A[(N+1):(2*N), 1:M] <- A[1:N, (M+1):(2*M)]
  A[(N+1):(2*N),(M+1):(2*M)] <- -D2*ydiff^2 - D1
  return(A)
}
#' @export
evaluate_DFtan_DFpart <- function(coef, Xe, X, Xb, Nb, kerngen){
  A <- DF_rbf2D(Xe, X, kerngen)

  B <- DF_rbf2D(Xe, Xb, kerngen)
  m <- nrow(Xb)
  Dn <- rbind(Matrix::diag(Nb[,1]),
              Matrix::diag(Nb[,2]))

  B <- B %*% Dn
  v <- cbind(A,B) %*% coef
  n <- nrow(Xe)
  vx <- v[1:n]
  vy <- v[(n+1):length(v)]
  return(cbind(vx, vy))
}
#' @export
evaluate_DFtan_CFpart <- function(coef, Xe, X, kerngen){
  A <- CF_rbf2D(Xe, X, kerngen)
  m <- nrow(X)
  v <- A %*% coef[1:(2*m)]
  n <- nrow(Xe)
  vx <- v[1:n]
  vy <- v[(n+1):length(v)]
  return(cbind(vx, vy))
}

#' Perform Hodge-Helmholtz Decomposition
#' @param dx delta x
#' @param dy delta y
#' @param Xi (matrix: (n_interior,2)) X interior nodes
#' @param Xb (matrix: (n_boundary, 2)) X boundary nodes
#' @param Nb (matrix: (n_boundary, 2)) Normal vectors at boundary nodes
#' @param Xe (matrix: (n, 2)) X domain points at which to evaluate projection
#' @param kerngen kernel generator function
#' @return List containing 
#' \itemize{
#' \item{coef}{projection coefficients}
#' \item{divfree}{divergence free projection at Xe}
#' \item{curlfree}{curlfree projection at Xe}
#' }
#' @export
hh_proj <- function(dx,dy, Xi, Xb, Nb, Xe, kerngen){
  d <- matrix(c(dx, dy, rep(0, nrow(Xb))), ncol = 1)
  dfpg <- build_DFprojGram(Xi, Xi, Xb, Nb, kerngen)
  G <- Matrix::bdiag(dfpg$A, dfpg$A)
  G <- rbind(cbind(G, dfpg$B),
             cbind(Matrix::t(dfpg$B), dfpg$C))
  coef <- solve(G, d)
  vs <- evaluate_DFtan_DFpart(coef, Xe, Xi, Xb, Nb, kerngen)
  ws <- evaluate_DFtan_CFpart(coef, Xe, Xi, kerngen)
  # hs <- evaluate_CFnorm_CFpart(coef, )
  return(list(coef = coef, divfree = vs, curlfree = ws))
}
