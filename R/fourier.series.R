#' @title Truncated Fourier Series
#' @description Truncated Fourier Series representation for an even periodic function. It does not give appropriate results for odd functions or non-odd and non-even functions.
#'
#' @param f_t The periodic function to be expressed as a truncated fourier series
#' @param k Order of trucncation
#'
#' @return List containing the truncated fourier series and the fourier coefficients
#' @export
#'
#' @examples fourier.series(y_simulated, k = 2)
#' @examples fourier.series(y_simulated, k = 1)
fourier.series <- function (f_t, k = 1){

  # Series in which we check the results
  T = length(f_t)
  frequency = c(0,2*pi)
  freqx = seq(frequency[1], frequency[2], length.out = T)

  # Fourier coefficients
  A_k = rep(0,k)
  delta = freqx[2]-freqx[1]

  for(j in 1:k){
    A_k[j] = (2 / (frequency[2] - frequency[1])) * sum(f_t * cos(j * freqx) * delta) # Integral approximation
  }

  #Fourier Series
  A_0 = mean(f_t)
  f_t_fourier = rep(A_0,T)
  for(j in 1:k){
    freq_t = j*freqx
    f_t_fourier = f_t_fourier + A_k[j] * cos(freq_t)
  }

  answer = list(f_t_fourier, A_k)
  names(answer) <- c("tfs", "coef")

  return(answer)
}
