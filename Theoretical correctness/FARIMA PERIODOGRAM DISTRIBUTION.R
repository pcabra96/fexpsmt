
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
T = 2^15

for (sim in 1:NUMBER_OF_SIMULATIONS){

    ################################################################################
    # PARAMETERS OF THE SIMULATION
    ################################################################################
    ar_coef = runif(1, min = -0.9, max = 0.9)
    ma_coef = runif(1, min = -0.9, max = 0.9)
    if(abs(ar_coef-ma_coef)<0.05 | abs(ar_coef+ma_coef)<0.05){
      ma_coef = -runif(1, min = -0.9, max = 0.9)
    }
    d = runif(1, min = 0.01, max = 0.45)

    ################################################################################
    # SIMULATION
    ################################################################################
    y = sim.farima(ar = ar_coef, ma = ma_coef, d = d, T = T)

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

    spec = farima.spectrum(ar = ar_coef, ma = ma_coef, d = d, n.freq = T)

    I = yper/spec[1:(T/2-1)]

    fit1 <- fitdistr(I, "exponential")
    fitted_exp_rate[sim] = fit1$estimate
    goodess_of_fit[sim] = ks.test(I, "pexp", 1)[["p.value"]]
}

################################################################################
# VISUALIZE THE RESULTS
################################################################################
boxplot(fitted_exp_rate, ylim = c(0.9,1.1), main = "Boxplot of fitted exp(1) rate", ylab = "exp(1) rate")
abline(h=1, col = "red")

plot(goodess_of_fit, main = "Goodness of fit of I(w)/f(w) to an exp(1)", ylab = "p.value")
abline(h=0.01, col = "red")
sum(goodess_of_fit<0.01)/length(goodess_of_fit)
