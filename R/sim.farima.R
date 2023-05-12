#' Function to simulate a FARIMA(p,d,q) process from its theoretical spectral density.
#'
#' @param ar the autoregressive coefficients
#' @param ma the moving average coefficients
#' @param d the long memorry parameter
#' @param T the length of the time series to be simulated
#'
#' @return a time series of lenth T
#' @export
#'
#' @examples sim.farima(ar = 0.5)
#' @examples sim.farima(ar = c(0.5,0.2),T = 128)
#' @examples sim.farima(ar = 0.5, ma = -0.8, d =0.2, T = 8192)

sim.farima <- function(ar = 0, ma = 0, d = 0, T = 512){

  # Assign values to parameters
  i = complex(real = 0, imaginary = 1)
  ar_coef = ar
  ma_coef = ma
  d = d
  M = T*2
  frequency = c(0,pi)
  freq <- seq.int(frequency[1], frequency[2], length.out = M)
  spectral_ts = farima.spectrum(ar = ar_coef, ma = ma_coef, d = d, frequency = frequency, n.freq = T+1)

  # Approx f(0) when there is long memory component.
  if(d!=0){
    freq_2 = seq(frequency[1],frequency[2],length.out = T+1)
    lamda_aprox = (freq_2[2]-freq_2[1])/T
    H = d+0.5
    lamda_aprox = lamda_aprox^(1-2*H)
    c_f = (1/(2*pi))*sin(H*pi)*gamma(2*H+1)
    spectral_ts[1] = lamda_aprox*c_f
  }

  Ws = rnorm(M)
  Ws = matrix(Ws, ncol = 2)

  Vs = rep(0,M)

  # Simple
  Vs[1] = sqrt(spectral_ts[1])*Ws[1,1] # Variance
  Vs[(M/2+1)]= sqrt(spectral_ts[(M/2+1)])*Ws[1,2] # Maximum lag

  # Composed
  # Real part 1 <= j < M/2
  t_1 = seq(2,(M/2),1)
  Vs[t_1] = sqrt(0.5*spectral_ts[t_1])*(Ws[t_1,1]+i*Ws[t_1,2])

  # Complex conjugate
  Vs[2+M-t_1] = Conj(Vs[t_1])

  # Inverse fast fourier transform
  y_fft = sqrt(2*pi/M)*Re(fft(Vs, inverse = T))[1:(M/2)]

  return(y_fft)
}
