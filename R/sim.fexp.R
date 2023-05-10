sim.fexp <- function(ck, d = 0, T = 512){

  # Assign values to parameters
  i = complex(real = 0, imaginary = 1)
  d = d
  M = T*2
  frequency = c(0,2*pi)
  freq <- seq.int(frequency[1], frequency[2], length.out = M)

  spectral_ts = fexp.spectrum(ck = ck, d = 0, frequency = frequency, n.freq = M)

  # Include long memory component
  l.m <- sapply(freq, function(x) abs(2*sin(x*(1/2)))^(-2*d))
  spectral_ts = spectral_ts*l.m

  # Approx f(0) when there is long memory component.
  if(d!=0){
    freq_2 = seq(frequency[1],frequency[2],length.out = T+1)
    lamda_aprox = (freq_2[2]-freq_2[1])/T
    H = d+0.5
    lamda_aprox = lamda_aprox^(1-2*H)
    c_f = (1/(2*pi))*sin(H*pi)*gamma(2*H+1)
    spectral_ts[1] = lamda_aprox*c_f
  }
  spectral_ts = spectral_ts[1:(T+1)]

  # Assign values to parameters
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
  y_fft = Re(fft(Vs, inverse = T))[1:(M/2)]/sqrt(M)
  return(y_fft)
}
