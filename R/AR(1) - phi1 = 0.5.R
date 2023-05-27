################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

ar_coef = -0.9
T = 2^12
true_spectrum = farima.spectrum(ar = ar_coef, n.freq = T)
################################################################################
# SIMULATION
################################################################################

N_SIMULATIONS = 500
fit_own_coef = rep(0,N_SIMULATIONS)
fit_r_coef = rep(0,N_SIMULATIONS)
fit_own_exp = rep(0,N_SIMULATIONS)
fit_r_exp = rep(0,N_SIMULATIONS)
p_val_own_exp = rep(0,N_SIMULATIONS)
p_val_r_exp = rep(0,N_SIMULATIONS)

for (sim in 1:N_SIMULATIONS) {
  # AR OWN simulations
  y_own = sim.farima(ar = ar_coef, T = T)

  # AR EXISTING PACKAGE
  y_r = arima.sim(model = list(ar = ar_coef), n = T)

  ################################################################################
  # FITTING
  ################################################################################

  fit_own_coef[sim] = fit.farima(y_own, p = 1)[["ar"]]
  fit_r_coef[sim] = fit.farima(y_r, p = 1)[["ar"]]

  ################################################################################
  # PERIODOGRAM
  ################################################################################

  n = length(y_own)

  # Fundamental frequencies
  mhalfm <- (n-1) %/% 2L
  w <- 2*pi/n * (1:mhalfm)

  # Periodogram ordinates by FFT of own simulation
  per_own = (Mod(fft(y_own))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper_own = per_own[2: ((n+1) %/% 2)]
  I_own = yper_own/true_spectrum[1:(T/2-1)]

  # Periodogram ordinates by FFT of R simulation
  per_r = (Mod(fft(y_r))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper_r = per_r[2: ((n+1) %/% 2)]
  I_r = yper_r/true_spectrum[1:(T/2-1)]

  # Periodogram distribution of own simulation
  fit_own_exp[sim] <- fitdistr(I_own, "exponential")[["estimate"]]
  p_val_own_exp[sim] = ks.test(I_own, "pexp", 1)[["p.value"]]

  # Periodogram distribution of R simulation
  fit_r_exp[sim] <- fitdistr(I_r, "exponential")[["estimate"]]
  p_val_r_exp[sim] = ks.test(I_r, "pexp", 1)[["p.value"]]
}

# Fitted AR coefficient
par(mfrow=c(1,2))
boxplot(fit_own_coef, main = "Fitted AR(1) own simulation")
abline(h=ar_coef, col = "red")
boxplot(fit_own_coef, main = "Fitted AR(1) R simulation")
abline(h=ar_coef, col = "red")

# Fitted AR periodogram
par(mfrow=c(1,2))
boxplot(fit_own_exp, main = "Fitted EXP(1) own simulation")
abline(h=1, col = "red")
boxplot(fit_r_exp, main = "Fitted EXP(1) R simulation")
abline(h=1, col = "red")

# P value
par(mfrow=c(1,2))
plot(p_val_own_exp, main = "p.value own simulation")
abline(h=0.05, col = "red")
plot(p_val_r_exp, main = "p.value R simulation")
abline(h=0.05, col = "red")

# Check figures of last simulation
par(mfrow=c(2,1))
plot(y_own, main = "Own simulation", type = "l")
plot(y_r, main = "R simulation", type = "l")

par(mfrow=c(1,1))
plot(yper_r, main = "Own periodogram")
lines(true_spectrum, type = "l", col = "red")
plot(yper_r, main = "R periodogram")
lines(true_spectrum, type = "l", col = "red")







