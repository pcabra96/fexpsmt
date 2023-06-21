################################################################################
##----------------------------------------------------------------------------##
## INDEX                                                                      ##
##----------------------------------------------------------------------------##
################################################################################

# 1. PACKAGES
# 2. SEED
# 3. SIMULATION PARAMETERS
# 4. SIMULATION
# 5. RESULTS
# 5.1. TIME
# 5.2. TIME DOMAIN PARAMETER
# 5.3. FREQUENCY DOMAIN PARAMETER
# 5.4. FREQUENCY DOMAIN GOODNESS OF FIT

N_SIMULATIONS = 1000
T = 2^13
ar_coef = 0.7
K = c(1,2,3,4,5,6,7,8)

fit_own_coef = matrix(0,nrow = N_SIMULATIONS, ncol = length(K))
fit_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(K))
p_val_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(K))

start = Sys.time()

for (sim in 1:N_SIMULATIONS) {
  for (j in 1:length(K)) {

  true_spectrum = farima.spectrum(ar = ar_coef, n.freq = T)
  fourier_series = fourier.series(f_t = log(true_spectrum), k = K[j])
  coef = fourier_series$coef

  ################################################################################
  # SIMULATION
  ################################################################################

  y = sim.fexp(ck = coef, T = T)

  ################################################################################
  # FITTED VALUE
  ################################################################################

  fit_own_coef[sim,j] = fit.farima(y, p = 1)[["ar"]]

  ################################################################################
  # PERIODOGRAM
  ################################################################################

  n = length(y)

  # Fundamental frequencies
  mhalfm <- (n-1) %/% 2L
  w <- 2*pi/n * (1:mhalfm)

  # Periodogram ordinates by FFT of own simulation
  per_own = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper_own = per_own[2: ((n+1) %/% 2)]
  I_own = yper_own/true_spectrum[1:(T/2-1)]

  # Periodogram distribution of own simulation
  fit_own_exp[sim,j] = fitdistr(I_own, "exponential")[["estimate"]]
  p_val_own_exp[sim,j] = ks.test(I_own, "pexp", 1)[["p.value"]]
  }
}

end = Sys.time()
print(end-start)

boxplot(fit_own_coef, ylim = c(ar_coef-0.2,ar_coef+0.2))
abline(h = ar_coef, col = "red")

boxplot(fit_own_exp, ylim = c(1-0.2,1+0.2))
abline(h = 1, col = "red")

boxplot(p_val_own_exp)
abline(h = 0.05, col = "red")
