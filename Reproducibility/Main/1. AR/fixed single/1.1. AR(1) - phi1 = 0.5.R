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

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(forecast)
require(MASS)
library(latex2exp)

################################################################################
##----------------------------------------------------------------------------##
## 2. SEED                                                                    ##
##----------------------------------------------------------------------------##
################################################################################

set.seed(1)

################################################################################
##----------------------------------------------------------------------------##
## 3. SIMULATION PARAMETERS                                                   ##
##----------------------------------------------------------------------------##
################################################################################

PROCESS = "AR"
symbol = "\\phi"
coef = 0.5
POWER = 7:14
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION                                                              ##
##----------------------------------------------------------------------------##
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
  true_spectrum = farima.spectrum(ar = coef, n.freq = T)
  for (sim in 1:N_SIMULATIONS) {

    # AR OWN simulations
    own_start = Sys.time()
    y_own = sim.farima(ar = coef, T = T)
    own_end = Sys.time()
    time_own[sim,j] = own_end - own_start

    # AR EXISTING PACKAGE
    r_start = Sys.time()
    y_r = arima.sim(model = list(ar = coef), n = T)
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
##----------------------------------------------------------------------------##
## 5. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/1. ",PROCESS,"/")

################################################################################
# 5.1. TIME
################################################################################

#--------------------------------------------
# 5.1.1 TIME BOXPLOTS
#--------------------------------------------

# Limits
lim_inf = min(c(time_own,time_r))
lim_sup = max(c(time_own,time_r))

#Boxplots
par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
boxplot(time_own, ylim=c(0,0.02), names = names, ylab = "time (s)", xlab = "T")
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
boxplot(time_r,ylim=c(0,0.02), names = names, ylab = "time (s)", xlab = "T")
title(main = "stats", cex.main = 0.8, line = 0.5)
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 1.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

#--------------------------------------------
# 5.1.2 TIME AVERAGES
#--------------------------------------------

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = colMeans(time_own), col = "blue", type = "o", ylim=c(lim_inf,lim_sup/2), ylab = "time (s)", labels = FALSE, xlab = "T")
lines(x = POWER, y = colMeans(time_r), col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)
main = paste0("Average running time for $\\",N_SIMULATIONS,"\\{ y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 2.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# 5.2. TIME DOMAIN PARAMETER
################################################################################

#--------------------------------------------
# 5.2.1 TIME DOMAIN COEFFICIENTS BOXPLOTS
#--------------------------------------------

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{",symbol,"_1}$")
boxplot(fit_own_coef, names = names, xlab = "T", ylab=TeX(main))
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
abline(h=coef, col = "red")
main = paste0("$\\hat{",symbol,"_1}$")
boxplot(fit_r_coef, names = names, xlab = "T", ylab = TeX(main))
title(main = "stats", cex.main = 0.8, line = 0.5)
abline(h=coef, col = "red")
main = paste0("Fitted $",symbol,"_1$ for$\\ ",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 3.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)

#--------------------------------------------
# 5.2.2 COEFFICIENTS AVERAGES AND DISPERSION
#--------------------------------------------

# Calculate standard deviation for each column
mean_phi_mse_own <- colMeans((fit_own_coef - coef)^2)
mean_phi_mse_r <- colMeans((fit_r_coef - coef)^2)
sd_phi_mse_own <- apply(((fit_own_coef - coef)^2), 2, sd)
sd_phi_mse_r <- apply(((fit_r_coef - coef)^2), 2, sd)

# Plot mean with error bars for 'fit_own_phi'
par(mfrow = c(1, 1), mar=c(5,5,3,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{",symbol,"}_1$ MSE")
plot(x = POWER, y = mean_phi_mse_own, col = "blue", type = "o", ylab = TeX(main), labels = FALSE, xlab = "T", ylim = c(min(mean_phi_mse_own-sd_phi_mse_own),max(mean_phi_mse_own+sd_phi_mse_own)))
lines(x = POWER, y = mean_phi_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col = "black")
legend("topright", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)

for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i], y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] - sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] + sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i], y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] - sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] + sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
}
main = paste0("MSE of $\\hat{",symbol,"}_1$ with $\\",N_SIMULATIONS," \\{y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 4.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# 5.3. FREQUENCY DOMAIN PARAMETER
################################################################################

#--------------------------------------------
# 5.3.1 FREQUENCY DOMAIN COEFFICIENTS BOXPLOTS
#--------------------------------------------

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{lambda}_{MLE}$")
boxplot(fit_own_exp, names = names, xlab = "T", ylab = TeX(main))
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
abline(h=1, col = "red")
main = paste0("$\\hat{lambda}_{MLE}$")
boxplot(fit_r_exp, names = names, xlab = "T", ylab =TeX(main))
title(main = "stats", cex.main = 0.8, line = 0.5)
abline(h=1, col = "red")
main = paste0("Fitted $\\lambda$ for ", N_SIMULATIONS, " $\\{I(\\omega_k)^*_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{k=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 5.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)

#--------------------------------------------
# 5.3.2 FREQUENCY DOMAIN COEFFICIENTS BOXPLOTS
#--------------------------------------------

mean_own_exp <- colMeans((fit_own_exp-1)^2)
mean_r_exp <- colMeans((fit_r_exp-1)^2)
sd_own_exp <- apply((fit_own_exp-1)^2, 2, sd)
sd_r_exp <- apply((fit_r_exp-1)^2, 2, sd)

min_lim = min(c((mean_own_exp-sd_own_exp),(mean_r_exp-sd_r_exp)))
max_lim = max(c((mean_own_exp+sd_own_exp),(mean_r_exp+sd_r_exp)))

# Plot means with error bars
par(mfrow=c(1,1), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_own_exp, col = "blue", type = "o", ylab = TeX("MSE $\\hat{lambda}$"), xlab = "T", ylim = c(min_lim, max_lim), labels = FALSE)
lines(x = POWER, y = mean_r_exp, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)

legend("topright", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)
abline(h=1, col = "black")

# Add error bars
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_own_exp[i] - sd_own_exp[i], x1 = POWER[i], y1 = mean_own_exp[i] + sd_own_exp[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp[i] - sd_own_exp[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp[i] - sd_own_exp[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp[i] + sd_own_exp[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp[i] + sd_own_exp[i], col = "blue")
  segments(x0 = POWER[i], y0 = mean_r_exp[i] - sd_r_exp[i], x1 = POWER[i], y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_r_exp[i] - sd_r_exp[i], x1 = POWER[i] + 0.1, y1 = mean_r_exp[i] - sd_r_exp[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_r_exp[i] + sd_r_exp[i], x1 = POWER[i] + 0.1, y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
}

main = paste0("MSE of fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{",PROCESS,"(",symbol,"_1=",coef,")_t}(\\omega_k)\\}_{k=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)
graph_name = "Figure 6.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# 5.4. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

#--------------------------------------------
# 5.4.1 GOODNESS OF FIT BOXPLOTS
#--------------------------------------------

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
boxplot(p_val_own_exp, ylab = "p.value", names = names, xlab = "T")
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
abline(h=0.05, col = "red")
boxplot(p_val_r_exp, ylab = "p.value", names = names, xlab = "T")
title(main = "stats", cex.main = 0.8, line = 0.5)
abline(h=0.05, col = "red")
main = paste0("$H_0: \\ \\{I(\\omega_k)^*_{",PROCESS,"(",symbol,"_1 = ", coef, ")_t}\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 7.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# Last checks
################################################################################

# Check figures of last simulation
par(mfrow=c(2,1))
main = paste0("One realization of $\\{y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{",2,"^{",POWER[length(POWER)],"}}$")
plot(y_own, main = TeX(main), type = "l", ylab = TeX("$y_{fexpmst}$"), xlab = "")
main = paste0("One realization of $\\{y_{",PROCESS,"(",symbol,"_1=",coef,")_t}\\}_{t=1}^{",2,"^{",POWER[length(POWER)],"}}$")
plot(y_r, main = TeX(main), type = "l", ylab = TeX("$y_{stats}$"), xlab = "")

graph_name = "Figure 8.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

# par(mar = c(bottom, left, top, right)
par(mfrow=c(2,1), mar = c(3, 5, 2, 2))
main = paste0("One realization of $\\{I^*_{",PROCESS,"(",symbol,"_1 = ", coef, ")_t}(\\omega_k)\\}_{k=1}^{",2,"^{",POWER[length(POWER)],"}-1}$")
plot(yper_own, ylab = TeX("$I(\\omega)_{fexpmst}$"), main = TeX(main), ylim = c(0,6))
lines(true_spectrum, type = "l", col = "red")
main = paste0("One realization of $\\{I(\\omega_k)^*_{",PROCESS,"(",symbol,"_1 = ", coef, ")_t}(\\omega_k)\\}_{k=1}^{",2,"^{",POWER[length(POWER)],"}-1}$")
plot(yper_r, ylab = TeX("$I(\\omega)_{stats}$"), main = TeX(main), ylim = c(0,6))
lines(true_spectrum, type = "l", col = "red")

graph_name = "Figure 9.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()
