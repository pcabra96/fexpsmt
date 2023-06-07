#' Fit the FARIMA(p,d,q) parameters using Whittle estimation
#'
#' @param y Time Series
#' @param max.p Integer representing the order of the autoregressive component of the process
#' @param max.q Integer representing the order of the moving average component of the process
#' @param max.d 1 if the process has long memory and 0 if it does not
#'
#' @return a list containing the estimated parameters of a FARIMA(p,d,q) process
#' @export
#'
#' @examples fit.farima(y, max.p = 1)
#' @examples fit.farima(y, max.p = 1, max.q = 1, max.d = 1)

fit.farima <- function(y, p = 0, q = 0, d = 0){
  n = length(y)
  # Fundamental frequencies
  mhalfm <- (n-1) %/% 2L
  w <- 2*pi/n * (1:mhalfm)

  # Periodogram ordinates by FFT
  per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper = per[2: ((n+1) %/% 2)]

  # Optimization parameters
  p = p
  q = q
  d_long = d
  x = w
  start_vals = c(runif(n = p+q, min = -0.5, max = 0.5),0)

  obj.Whittle = function(theta, p = 0, q = 0, d_long=0, x){

    I_w = yper
    f_w = spectral.density.farima(theta, p = p, q = q, d_long = d_long, x = x)

    SPO = I_w/f_w

    return(sum(SPO))
  }

  # Define the lower and upper bounds for the optimization
  lower = c(rep(-0.9999, p+q),0)
  upper = c(rep(0.9999, p+q),0.55)

  opt_res = optim(par = start_vals, fn = obj.Whittle, p = p, q = q, d_long = d_long, x = x, lower = lower, upper = upper, method = "L-BFGS-B")
  # Print the optimized parameters
  if(p>0) {
    ar_parameters = opt_res$par[1:p]
  }
  else{
    ar_parameters = 0
  }
  if(q>0){
    ma_parameters = opt_res$par[(p+1):(p+q)]
  }
  else{
    ma_parameters = 0
  }
  if(d_long>0){
    d_parameter = opt_res$par[p+q+1]
  }
  else{
    d_parameter = 0
  }
  answer = list(ar_parameters, ma_parameters, d_parameter)
  names(answer) <- c("ar", "ma", "d")

  return(answer)
}
