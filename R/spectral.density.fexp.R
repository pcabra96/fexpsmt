#' Title
#'
#' @param eta
#' @param p
#' @param d_long
#' @param x
#'
#' @return
#' @export
#'
#' @examples

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
