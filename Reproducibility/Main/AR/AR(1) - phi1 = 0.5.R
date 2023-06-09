################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)

################################################################################
# SEED
################################################################################

set.seed(1)

################################################################################
# PARAMETERS
################################################################################

ar_coef = 0.5
POWER = 7:14
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# SIMULATION
################################################################################

N_SIMULATIONS = 1000
fit_own_coef = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_coef = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
time_own = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
time_r = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

for (j in 1:length(POWER)) {
  T = 2^(POWER[j])
  true_spectrum = farima.spectrum(ar = ar_coef, n.freq = T)
  for (sim in 1:N_SIMULATIONS) {

    # AR OWN simulations
    own_start = Sys.time()
    y_own = sim.farima(ar = ar_coef, T = T)
    own_end = Sys.time()
    time_own[sim,j] = own_end - own_start

    # AR EXISTING PACKAGE
    r_start = Sys.time()
    y_r = arima.sim(model = list(ar = ar_coef), n = T)
    r_end = Sys.time()
    time_r[sim,j] = r_end - r_start

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

# Run time AR coefficient
par(mfrow=c(1,2), mar = c(5, 4, 2, 2)) # mar = c(bottom, left, top, right))
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{AR(phi_1=",ar_coef,")_t,own}\\}_{t=1}^{T}$")
boxplot(time_own, ylim=c(0,0.02), main = TeX(main), names = names, ylab = "time (s)", xlab = "T")
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{AR(phi_1=",ar_coef,")_t,R}\\}_{t=1}^{T}$")
boxplot(time_r,ylim=c(0,0.02), main = TeX(main), names = names, ylab = "time (s)", xlab = "T")
graph_name = "Figure 1.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)

# Fitted AR coefficient
par(mfrow=c(1,2), mar = c(4.5, 5, 3, 1)) # mar = c(bottom, left, top, right))
main = paste0("Fitted $\\phi_1$ for$\\ ",N_SIMULATIONS," \\ \\{y_{AR(phi_1=",ar_coef,")_t,own}\\}_{t=1}^{T}$")
boxplot(fit_own_coef, main = TeX(main), names = names, xlab = "T", ylab=TeX("$\\hat{phi_1}$"))
abline(h=ar_coef, col = "red")
main = paste0("Fitted $\\phi_1$ for$\\ ",N_SIMULATIONS," \\ \\{y_{AR(phi_1=",ar_coef,")_t,R}\\}_{t=1}^{T}$")
boxplot(fit_r_coef, main = TeX(main), names = names, xlab = "T", ylab = TeX("$\\hat{phi_1}$"))
abline(h=ar_coef, col = "red")
graph_name = "Figure 2.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)

# Fitted AR periodogram
par(mfrow=c(1,2), mar = c(4.5, 5, 2, 2))
main = paste0("Fitted $\\lambda$ for ", N_SIMULATIONS, " $\\{I(\\omega_k)^*_{AR(phi_1=",ar_coef,")_t,own}\\}_{k=1}^{T-1}$")
boxplot(fit_own_exp, main = TeX(main), names = names, xlab = "T", ylab = TeX("$\\hat{lambda}_{MLE}$"))
abline(h=1, col = "red")
main = paste0("Fitted $\\lambda$ for ", N_SIMULATIONS, " $\\{I(\\omega_k)^*_{AR(phi_1=",ar_coef,")_t,R}\\}_{k=1}^{T-1}$")
boxplot(fit_r_exp, main = TeX(main), names = names, xlab = "T", ylab =TeX("$\\hat{lambda}_{MLE}$"))
abline(h=1, col = "red")
graph_name = "Figure 3.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)

# P value
par(mfrow=c(1,2), mar = c(4.5, 5, 2, 2))
main = paste0(N_SIMULATIONS," $H_0: \\ \\{I(\\omega_k)^*_{AR(phi_1 = ", ar_coef, ")_t,own}\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
boxplot(p_val_own_exp, main = TeX(main), ylab = "p.value", names = names, xlab = "T")
abline(h=0.05, col = "red")
main = paste0(N_SIMULATIONS," $H_0: \\ \\{I(\\omega_k)^*_{AR(phi_1 = ", ar_coef, ")_t,R}\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
boxplot(p_val_r_exp, main = TeX(main), ylab = "p.value", names = names, xlab = "T")
abline(h=0.05, col = "red")
graph_name = "Figure 4.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)

# Check figures of last simulation
par(mfrow=c(2,1))
main = paste0("One realization of $\\{y_{AR(phi_1=",ar_coef,")_t,own}\\}_{t=1}^{",2,"^{",POWER[length(POWER)],"}}$")
plot(y_own, main = TeX(main), type = "l", ylab = TeX("$y_{own}$"), xlab = "")
main = paste0("One realization of $\\{y_{AR(phi_1=",ar_coef,")_t,R}\\}_{t=1}^{",2,"^{",POWER[length(POWER)],"}}$")
plot(y_r, main = TeX(main), type = "l", ylab = TeX("$y_{R}$"), xlab = "")
graph_name = "Figure 5.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)

# par(mar = c(bottom, left, top, right)
par(mfrow=c(2,1), mar = c(3, 5, 2, 2))
main = paste0("One realization of $\\{I(\\omega_k)^*_{AR(\\phi_1 = ", ar_coef, ")_t,own}\\}_{k=1}^{",2,"^{",POWER[length(POWER)],"}-1}$")
plot(yper_own, ylab = TeX("$I(\\omega)$"), main = TeX(main), ylim = c(0,6))
lines(true_spectrum, type = "l", col = "red")
main = paste0("One realization of $\\{I(\\omega_k)^*_{AR(\\phi_1 = ", ar_coef, ")_t,R}\\}_{k=1}^{",2,"^{",POWER[length(POWER)],"}-1}$")
plot(yper_r, ylab = TeX("$I(\\omega)$"), main = TeX(main), ylim = c(0,6))
lines(true_spectrum, type = "l", col = "red")
graph_name = "Figure 6.png"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/", graph_name)
dev.print(device = png, filename = path, width = 1800, height = 1100, res=200)
