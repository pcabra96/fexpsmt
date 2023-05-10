#' @title Spectral density of an arbitrary FEXP process.
#' @description It is used to fit the simulated time series.
#' @param eta Vector with the parameters
#' @param p Integer representing the order of the FEXP process
#' @param d_long Integer with zero or one representing if one considers the long memory component
#' @param x Frequencies to be used
#'
#' @return Spectral density of the arbitrary FARIMA process
#' @export
#'
#' @examples spectral.density.farima(eta = c(0.5,0.5,0.3), p = 2, d_long=1, x)

spectral.density.fexp <- function(eta = c(), p = 0, d_long=0, x){

  if(p > 0) {
    ck <- eta[1:(p)]
    px = outer(x, 1:p)
    Rar = cos(px) %*% ck

  } else {
    phi <- numeric(0)
    far <- 1
  }

  if(d_long > 0){
    d <- eta[p+1]
    lm = sapply(x, function(y) abs(2*sin(y*(1/2)))^(-2*d))
  }
  else{
    psi = 0
    lm <- 1
  }
  spec <- lm*exp(Rar)

  spectrum = (1/(2*pi))*spec

  return(spectrum)
}
