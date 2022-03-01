#' Get Central Differences
#' @description Calculate central differences (e.g. to approximate a first 
#' derivative of f(x))
#' @param f function to cacluate differences of
#' @param x point in domain of function at which to calculate central difference
#' @param del delta of central difference (f(x+0.5del) - f(x-0.5del))/2
#'
#' @export
central_diff <- function(f, x, del){
  (f(x + 0.5*del) - f(x - 0.5*del))/del
}


#' Calculate gradient of f(x,y)
#' @param xy (matrix (n,2)) of points at which to calculate gradient of f
#'
#' @param f (function) scalar function with two real-valued inputs
#' @param del (scalar) delta to use in central difference approximation of 
#' derivatives (see \code{central_diff()})
#'
#' @export
calc_grad <- function(xy, f, del = 1e-5){
  vx <- purrr::map_dbl(1:nrow(xy), ~ central_diff(function(x)
    f(x, xy[.x, 2]), xy[.x, 1], del = del))
  vy <- purrr::map_dbl(1:nrow(xy), ~ central_diff(function(y)
    f(xy[.x, 1], y), xy[.x, 2], del = del))
  return(cbind(vx, vy))
}

#' Gradient of Gaussian density
#' @inheritParams central_diff
#' @param mu (2-dim vector) bi-variate normal mean
#' @param sigma (matrix(2,2)) bi-variate normal covariance
#' @export
gauss_grad <- function(xy, mu, sigma, del = 1e-5){
  vxvy <- calc_grad(xy, function(x,y){
    mvtnorm::dmvnorm(c(x,y), mu, sigma)
  }, del = del)
}
