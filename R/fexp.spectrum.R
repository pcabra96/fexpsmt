#' @title FEXP(p,d) theoretical spectral density
#' @description
#'
#' @param ck The p coefficients corresponding to the exp component of the model. ATTENTION: this ck is 2 times the coefficient proposed by Bloomfield in his original paper.
#' @param d The long memory component, that should be between 0 and 0.5
#' @param var.noise Variance of the process
#' @param n.freq Length of the simulation
#' @param frequency # An interval of the domain in which the theoretical spectral density is evaluated
#'
#' @return Theoretical spectral density of a FEXP process
#' @export
#'
#' @examples fexp.spectrum(ck = 1)
#' @examples fexp.spectrum(ck = c(2,1), n.freq = 2048)

fexp.spectrum <- function (ck = 0, d = 0, var.noise = 1, n.freq = 512, frequency=c(0,2*pi)){

  fexp_order <- length(ck)
  frequency = frequency
  freq <- seq.int(frequency[1], frequency[2], length.out = n.freq)

  if(fexp_order > 0) {
    ck = ck
    px = outer(freq, 1:fexp_order)
    Rar = cos(px) %*% ck
  }
  else {
    ck = numeric(0)
    far = 1
  }

  if(d>0){
    d = d
    lm = sapply(freq, function(y) abs(2*sin(y*(1/2)))^(-2*d))
  }
  else{
    psi = 0
    lm = 1
  }

  spec = lm * exp(Rar)
  spectrum = (1/(2*pi))*spec

  return(spectrum)
}
