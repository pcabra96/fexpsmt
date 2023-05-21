
################################################################################
# PACKAGES TO CHECK SIMULATIONS PROPERTIES
################################################################################

require(MASS)

################################################################################
# METAPARAMETERS
################################################################################
set.seed(0)
NUMBER_OF_SIMULATIONS = 1000
fitted_exp_rate = rep(0,NUMBER_OF_SIMULATIONS)
goodess_of_fit = rep(0,NUMBER_OF_SIMULATIONS)
T = 2^16

for (sim in 1:NUMBER_OF_SIMULATIONS){

  ################################################################################
  # PARAMETERS OF THE SIMULATION
  ################################################################################
  c_k = runif(1,min = -6, max = 6)
  d = 0

  ################################################################################
  # SIMULATION
  ################################################################################
  y = sim.fexp(ck = c_k, d = d, T = T)

  ################################################################################
  # PERIODOGRAM OF SIMULATED TIME SERIES
  ################################################################################
  n = length(y)

  # Fundamental frequencies
  mhalfm <- (n-1) %/% 2L
  w <- 2*pi/n * (1:mhalfm)

  # Periodogram ordinates by FFT
  per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
  yper = per[2: ((n+1) %/% 2)]

  spec = fexp.spectrum(ck = c_k, d = d, n.freq = T)

  I = yper/spec[1:(T/2-1)]

  fit1 <- fitdistr(I, "exponential")
  fitted_exp_rate[sim] = fit1$estimate
  goodess_of_fit[sim] = ks.test(I, "pexp", 1)[["p.value"]]
}

################################################################################
# VISUALIZE THE RESULTS
################################################################################

boxplot(fitted_exp_rate, ylim = c(0,1.5), main = "Boxplot of fitted exp(1) rate", ylab = "exp(1) rate")
abline(h=1, col = "red")

plot(goodess_of_fit, main = "Goodness of fit of I(w)/f(w) to an exp(1)", ylab = "p.value")
abline(h=0.01, col = "red")
sum(goodess_of_fit<0.01)/length(goodess_of_fit)


################################################################################
# Fiting an EXP(1)
################################################################################

vec = seq(-10,10,0.005)
n = length(vec)
coef = rep(0,n)

for (i in 1:n) {
  y = sim.fexp(ck = vec[i], T = 2^13)
  coef[i] = fit.fexp(y, p = 1, d = 0)[["c_k"]]
}
plot(x = vec, y = coef, type = "l", xlab = "true ck", ylab = "ck_hat", main = "True ck vs ckhat")
lines(x = vec, y = vec, type = "l", col = "red")

plot((coef-vec)^2)
abline
vec[3000]

################################################################################
y = sim.fexp(ck = c(-4,-2), T = 2^13)
plot(y, type = "l")
n = length(y)

# Fundamental frequencies
mhalfm <- (n-1) %/% 2L
w <- 2*pi/n * (1:mhalfm)
per = (Mod(fft(y))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
yper = per[2: ((n+1) %/% 2)]

Y = as.matrix(log(yper))
X = as.matrix(cbind(1, cos(w), cos(2*w)))

BETAS = solve(t(X) %*% X) %*% t(X) %*% Y
BETAS[c(2,3)]

