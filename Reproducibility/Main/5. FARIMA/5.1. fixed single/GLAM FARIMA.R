################################################################################
# PACKAGES
################################################################################

library(arfima)
library(fracdiff)
library(forecast)
require(MASS)
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

T = 2^14
ar_coef = 0.5
ma_coef = -0.2
d_coef  = 0.3
true_spectrum = farima.spectrum(ar = ar_coef,ma= ma_coef, d = d_coef, n.freq = T)
true_spectrum[1] = true_spectrum[length(true_spectrum)]

################################################################################
# OWN SIMULATION
################################################################################

N_SIMULATIONS = 1000

sigma_fexpmst_vec = rep(0,N_SIMULATIONS)
d_fexpmst_vec = rep(0,N_SIMULATIONS)
ar_coef_fexpmst_vec = rep(0,N_SIMULATIONS)
ma_coef_fexpmst_vec = rep(0,N_SIMULATIONS)

sigma_fracdiff_vec = rep(0,N_SIMULATIONS)
d_fracdiff_vec = rep(0,N_SIMULATIONS)
ar_coef_fracdiff_vec = rep(0,N_SIMULATIONS)
ma_coef_fracdiff_vec = rep(0,N_SIMULATIONS)

sigma_arfima_vec = rep(0,N_SIMULATIONS)
d_arfima_vec = rep(0,N_SIMULATIONS)
ar_coef_arfima_vec = rep(0,N_SIMULATIONS)
ma_coef_arfima_vec = rep(0,N_SIMULATIONS)

for (sim in 1:N_SIMULATIONS) {
    # 1. OWN SIMULATION
    y_own = sim.farima(ar = ar_coef, ma = ma_coef, d = d_coef, T = T)

    # 2. FRACDIFF
    y_fracdiff = fracdiff.sim(n = T, ar = ar_coef, ma = -ma_coef, d = d_coef)[["series"]]

    # 3. ARFIMA
    y_arfima = arfima.sim(n = T, list(phi = numeric(ar_coef), theta = numeric(ma_coef), dfrac = d_coef))

    ################################################################################
    # PERIODOGRAM
    ################################################################################

    n = length(y_own)

    # Fundamental frequencies
    mhalfm <- (n-1) %/% 2L
    w <- 2*pi/n * (1:mhalfm)

    # 1. OWN SIMULATION
    per_own = (Mod(fft(y_own))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
    yper_own = per_own[2: ((n+1) %/% 2)]
    I_own = yper_own/true_spectrum[1:(T/2-1)]

    # 2. FRACDIFF
    per_fracdiff = (Mod(fft(y_fracdiff))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
    yper_fracdiff = per_fracdiff[2: ((n+1) %/% 2)]
    I_fracdiff = yper_fracdiff/true_spectrum[1:(T/2-1)]

    # 3. ARFIMA
    per_arfima = (Mod(fft(y_arfima))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
    yper_arfima = per_arfima[2: ((n+1) %/% 2)]
    I_arfima = yper_arfima/true_spectrum[1:(T/2-1)]

    ################################################################################
    # GLM FITTING
    ################################################################################

    Y = as.matrix(cbind(log(yper_own), log(yper_fracdiff), log(yper_arfima)))
    X = as.matrix(cbind(1,log(abs(2*(sin(w/2)))),2*cos(w),cos(2*w)))

    BETA = solve(t(X) %*% X) %*% t(X) %*% Y

    # 1. OWN SIMULATION
    beta_0_fexpmst = BETA[1,1]
    beta_1_fexpmst = BETA[2,1]
    beta_2_fexpmst = BETA[3,1]
    beta_3_fexpmst = BETA[4,1]

    sigma_hat_fexpmst = exp(beta_0_fexpmst)*2*pi
    d_hat_fexpmst = -beta_1_fexpmst/2
    phi_hat_fexpmst = (beta_3_fexpmst + beta_2_fexpmst^2)/(2*beta_2_fexpmst)
    theta_hat_fexpmst = (beta_3_fexpmst-beta_2_fexpmst^2)/(2*beta_2_fexpmst)

    sigma_fexpmst_vec[sim] = sqrt(sigma_hat_fexpmst)
    d_fexpmst_vec[sim] = d_hat_fexpmst
    ar_coef_fexpmst_vec[sim] = phi_hat_fexpmst
    ma_coef_fexpmst_vec[sim] = theta_hat_fexpmst

    # 2. FRACDIFF
    beta_0_fracdiff = BETA[1,2]
    beta_1_fracdiff = BETA[2,2]
    beta_2_fracdiff = BETA[3,2]
    beta_3_fracdiff = BETA[4,2]

    sigma_hat_fracdiff = exp(beta_0_fracdiff)*2*pi
    d_hat_fracdiff = -beta_1_fracdiff/2
    phi_hat_fracdiff = (beta_3_fracdiff + beta_2_fracdiff^2)/(2*beta_2_fracdiff)
    theta_hat_fracdiff = (beta_3_fracdiff - beta_2_fracdiff^2)/(2*beta_2_fracdiff)

    sigma_fracdiff_vec[sim] = sqrt(sigma_hat_fracdiff)
    d_fracdiff_vec[sim] = d_hat_fracdiff
    ar_coef_fracdiff_vec[sim] = phi_hat_fracdiff
    ma_coef_fracdiff_vec[sim] = theta_hat_fracdiff

    # 3. ARFIMA
    beta_0_arfima = BETA[1,3]
    beta_1_arfima = BETA[2,3]
    beta_2_arfima = BETA[3,3]
    beta_3_arfima = BETA[4,3]

    sigma_hat_arfima = exp(beta_0_arfima)*2*pi
    d_hat_arfima = -beta_1_arfima/2
    phi_hat_arfima = (beta_3_arfima+beta_2_arfima^2)/(2*beta_2_arfima)
    theta_hat_arfima = (beta_3_arfima-beta_2_arfima^2)/(2*beta_2_arfima)

    sigma_arfima_vec[sim] = sqrt(sigma_hat_arfima)
    d_arfima_vec[sim] = d_hat_arfima
    ar_coef_arfima_vec[sim] = phi_hat_arfima
    ma_coef_arfima_vec[sim] = theta_hat_arfima
}

# SIGMA
plot(sigma_fexpmst_vec, type = "l", col = "blue", ylim=c(0,1), ylab = TeX("$\\sigma$"))
lines(sigma_fracdiff_vec, type = "l", col = "green")
lines(sigma_arfima_vec, type = "l", col = "red")
abline(h=1, col = "black")
legend("bottomright", legend = c("fepxmst", "fracdiff","farima"), col = c("blue", "green","red"), lty = 1)

# d
plot(d_fexpmst_vec, type = "l", col = "blue", ylim=c(d_coef-0.5,d_coef+0.5), ylab = TeX("$\\hat{d}$"))
lines(d_fracdiff_vec, type = "l", col = "green")
lines(d_arfima_vec, type = "l", col = "red")
abline(h=d_coef, col = "black")
legend("topright", legend = c("fepxmst", "fracdiff","farima"), col = c("blue", "green","red"), lty = 1)

# AR
plot(ar_coef_fexpmst_vec, type = "l", col = "blue", ylim=c(ar_coef-0.5,ar_coef+0.5), ylab = TeX("$\\hat{\\phi}_1$"))
lines(ar_coef_fracdiff_vec, type = "l", col = "green")
lines(ar_coef_arfima_vec, type = "l", col = "red")
abline(h=ar_coef, col = "black")
legend("topright", legend = c("fepxmst", "fracdiff","farima"), col = c("blue", "green","red"), lty = 1)

# MA
min_lim = min(ma_coef_fexpmst_vec,ma_coef_fracdiff_vec,ma_coef_arfima_vec)
max_lim = max(ma_coef_fexpmst_vec,ma_coef_fracdiff_vec,ma_coef_arfima_vec)
plot(ma_coef_fexpmst_vec, type = "l", col = "blue", ylim=c(ma_coef-0.5,ma_coef+0.5), ylab = TeX("$\\hat{\\theta}_1$"))
lines(ma_coef_fracdiff_vec, type = "l", col = "green")
lines(ma_coef_arfima_vec, type = "l", col = "red")
abline(h=ma_coef, col = "black")
legend("topright", legend = c("fepxmst", "fracdiff","farima"), col = c("blue", "green","red"), lty = 1)


