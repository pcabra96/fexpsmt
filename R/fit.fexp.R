#' Fit the FEXP(p,d) parameters using Whittle estimation
#'
#' @param y Time Series
#' @param p Integer representing the order of the exp component of the process
#' @param d 1 if the process has long memory and 0 if it does not
#'
#' @return a list containing the estimated parameters of a FEXP(p,d) process
#' @export
#'
#' @examples fit.fexp(y, max.p = 1)
#' @examples fit.fexp(y, max.p = 1, max.d = 1)


fit.fexp <- function(y, p = 0, d = 0, which_method = "Nelder-Mead"){
  n = length(y)
  # Fundamental frequencies
  mhalfm <- (n-1) %/% 2L
  w <- 2*pi/n * (1:mhalfm)

  # Periodogram ordinates by FFT
  per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper = per[2: ((n+1) %/% 2)]

  # Optimization parameters
  p = p
  d_long = d
  x = w
  start_vals = c(rep(0, p),0)

  obj.Whittle = function(theta, p = 0, d_long=0, x){

    I_w = yper
    f_w = spectral.density.fexp(theta, p = p, d_long = d_long, x = x)

    SPO = I_w/f_w

    return(sum(SPO))
  }

  # Define the lower and upper bounds for the optimization
  lower = c(rep(-Inf, p), 0)
  upper = c(rep(Inf, p), 0.55)

  opt_res = optim(par = start_vals, fn = obj.Whittle, p = p, d_long = d_long, x = x, lower = lower, upper = upper, method = which_method)
  # Print the optimized parameters
  if(p>0) {
    c_k_parameters = opt_res$par[1:p]
  }
  else{
    c_k_parameters = 0
  }
  if(d_long>0){
    d_parameter = opt_res$par[p+1]
  }
  else{
    d_parameter = 0
  }
  answer = list(c_k_parameters, d_parameter)
  names(answer) <- c("c_k", "d")

  return(answer)
}
