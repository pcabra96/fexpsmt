#' Theoretical Spectral density of a FARIMA(p,d,q) process.
#'
#' It is based on an existing function from stats but modified to include long memory.
#'
#' @param ar The coefficients of the autoregressive part of the process
#' @param ma The coefficients of the moving average part of the process
#' @param d # The long memory component, that should be between 0 and 0.5
#' @param var.noise # The variance of the process
#' @param n.freq # Length of the simulation
#' @param frequency # An interval of the domain in which the theoretical spectral density is evaluated
#'
#' @return Theoretical spectral density of a FARIMA process
#' @export
#'
#' @examples farima.spectrum(ar = 0.5)
#' @examples farima.spectrum(ar = c(0.5,-0.3), n.freq = 1024)
#' @examples farima.spectrum(ar = 0.5, ma = -0.2, d = 0.3, n.freq = 8192)

farima.spectrum <- function (ar = 0, ma = 0, d = 0, var.noise = 1, n.freq = 512, frequency=c(0,2*pi)){
  frequency = frequency
  ar.poly <- c(1, -ar)
  z.ar <- base::polyroot(ar.poly)
  check <- 0
  if (any(abs(z.ar) <= 1)) {
    cat("WARNING: Model Not Causal", "\n")
    check <- check + 1}

  ma.poly <- c(1, ma)
  z.ma <- base::polyroot(ma.poly)
  if (any(abs(z.ma) <= 1)) {
    cat("WARNING: Model Not Invertible", "\n")
    check <- check + 1
  }
  if (check > 0)
    stop("Try Again")
  ar.order <- length(ar)
  ma.order <- length(ma)
  for (i in 1:ar.order) {
    if ((ar[1] == 0 && ar.order == 1) || (ma[1] == 0 && ma.order == 1))
      break
    if (any(abs(z.ar[i] - z.ma[1:ma.order]) < 0.1)) {
      cat("WARNING: Parameter Redundancy", "\n")
      break}
  }
  if(abs(d)>0.5){
    stop("Error: d must be in the [-0.5,0.5] range")
  }

  freq <- seq.int(frequency[1], frequency[2], length.out = n.freq)
  cs.ar <- outer(freq, 1:ar.order, function(x, y) cos(x * y)) %*% ar
  sn.ar <- outer(freq, 1:ar.order, function(x, y) sin(x * y)) %*% ar
  cs.ma <- outer(freq, 1:ma.order, function(x, y) cos(x * y)) %*% -ma
  sn.ma <- outer(freq, 1:ma.order, function(x, y) sin(x * y)) %*% -ma
  l.m <- sapply(freq, function(x) abs(2*sin(x*(1/2)))^(-2*d))
  spec <- (1/(2*pi))*var.noise * l.m * ((1 - cs.ma)^2 + sn.ma^2)/((1 - cs.ar)^2 + sn.ar^2)
  return(spec)
}
