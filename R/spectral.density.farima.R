#' @title Spectral density of an arbitrary FARIMA process.
#' @description It is used to fit the simulated time series.
#' @param eta Vector with the parameters
#' @param p Integer representing the order of the auto regressive component
#' @param q Integer representing the order of the moving average component
#' @param d_long Integer with zero or one representing if one considers the long memory component
#' @param x Frequencies to be used
#'
#' @return Spectral density of the arbitrary FARIMA process
#' @export
#'
#' @examples spectral.density.farima(eta = c(0.5,0.5,0.3), p = 1, q = 1, d_long=1, x)
spectral.density.farima <- function(eta = c(), p = 0, q = 0, d_long=0, x){

  if(p > 0) {
    phi <- eta[1:p]
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% phi
    Iar <- sin(px) %*% phi

    far <- (1-Rar)^2 + Iar^2
  } else {
    phi <- numeric(0)
    far <- 1
  }

  if(q > 0) {
    psi <- eta[(p+1):(p+q)]
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi

    fma <- (1+Rma)^2 + Ima^2
  } else {
    psi <- numeric(0)
    fma <- 1
  }
  if(d_long>0){
    d <- eta[p+q+1]
    lm = sapply(x, function(y) abs(2*sin(y*(1/2)))^(-2*d))
  }
  else{
    psi = 0
    lm = 1
  }
  spec <- lm*(fma/far)

  spectrum = (1/(2*pi))*spec

  return(spectrum)
}
