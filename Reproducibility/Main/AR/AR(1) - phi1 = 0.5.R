################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

ar_coef = 0.5
POWER = 7:14
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# SIMULATION
################################################################################

N_SIMULATIONS = 500
fit_own_coef = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_coef = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

for (j in 1:length(POWER)) {
  T = 2^(POWER[j])
  true_spectrum = farima.spectrum(ar = ar_coef, n.freq = T)
  for (sim in 1:N_SIMULATIONS) {

    # AR OWN simulations
    y_own = sim.farima(ar = ar_coef, T = T)

    # AR EXISTING PACKAGE
    y_r = arima.sim(model = list(ar = ar_coef), n = T)

    ################################################################################
    # FITTING
    ################################################################################

    fit_own_coef[sim,j] = fit.farima(y_own, p = 1)[["ar"]]
    fit_r_coef[sim,j] = fit.farima(y_r, p = 1)[["ar"]]

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
    fit_own_exp[sim,j] <- fitdistr(I_own, "exponential")[["estimate"]]
    p_val_own_exp[sim,j] = ks.test(I_own, "pexp", 1)[["p.value"]]

    # Periodogram distribution of R simulation
    fit_r_exp[sim,j] <- fitdistr(I_r, "exponential")[["estimate"]]
    p_val_r_exp[sim,j] = ks.test(I_r, "pexp", 1)[["p.value"]]
  }
}

################################################################################
# VISUALIZE RESULTS
################################################################################

# Fitted AR coefficient
par(mfrow=c(1,2), mar = c(4.5, 5, 2, 2)) # mar = c(bottom, left, top, right))
boxplot(fit_own_coef, main = "Fitted AR(1) own simulation", names = names, xlab = "T", ylab=TeX("$\\hat{phi_1}$"))
abline(h=ar_coef, col = "red")
boxplot(fit_r_coef, main = "Fitted AR(1) R simulation", names = names, xlab = "T", ylab = TeX("$\\hat{phi_1}$"))
abline(h=ar_coef, col = "red")
graph_name = "Figure 1.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1597, height = 987, res=200)

# Fitted AR periodogram
par(mfrow=c(1,2), mar = c(4.5, 5, 2, 2))
boxplot(fit_own_exp, main = "Fitted EXP(1) own simulation", names = names, xlab = "T", ylab = TeX("$\\hat{lambda}$"))
abline(h=1, col = "red")
boxplot(fit_r_exp, main = "Fitted EXP(1) R simulation", names = names, xlab = "T", ylab =TeX("$\\hat{lambda}$"))
abline(h=1, col = "red")
graph_name = "Figure 2.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print( device = png, filename = path, width = 987, height = 610)

# P value
par(mfrow=c(1,2))
boxplot(p_val_own_exp, main = TeX("$H_0: \ \\{y_{own_t}\\}_{t=1}^{T} \\sim exp(\\lambda=1)$"),ylab = "p.value", names = names)
abline(h=0.05, col = "red")
boxplot(p_val_r_exp, main = TeX("$H_0: \ \\{y_{R_t}\\}_{t=1}^{T} \\sim exp(\\lambda=1)$"), ylab = "p.value", names = names)
abline(h=0.05, col = "red")
graph_name = "Figure 3.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print( device = png, filename = path, width = 987, height = 610)

# Check figures of last simulation
par(mfrow=c(2,1))
plot(y_own, main = "Own simulation", type = "l", ylab = TeX("$y_{own}$"), xlab = "")
plot(y_r, main = "R simulation", type = "l", ylab = TeX("$y_{R}$"), xlab = "")
graph_name = "Figure 4.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print( device = png, filename = path, width = 987, height = 610)


# par(mar = c(bottom, left, top, right)
par(mfrow=c(2,1), mar = c(3, 5, 2, 2))
plot(yper_own, main = TeX("Own $I(\\omega)$"), ylab = TeX("$I_{AR(1)}(\\omega)$"))
lines(true_spectrum, type = "l", col = "red")
plot(yper_r, main = TeX("R $I(\\omega)$"), ylab = TeX("$I_{AR(1)}(\\omega)$"))
lines(true_spectrum, type = "l", col = "red")
graph_name = "Figure 5.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print( device = png, filename = path, width = 987, height = 610)
