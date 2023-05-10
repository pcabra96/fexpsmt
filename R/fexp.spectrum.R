#' Title
#'
#' @param ck
#' @param d
#' @param var.noise
#' @param n.freq
#' @param frequency
#'
#' @return
#' @export
#'
#' @examples

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
    lm = sapply(x, function(y) abs(2*sin(y*(1/2)))^(-2*d))
  }
  else{
    psi = 0
    lm = 1
  }

  spec = lm * exp(Rar)
  spectrum = (1/(2*pi))*spec

  return(spectrum)
}
