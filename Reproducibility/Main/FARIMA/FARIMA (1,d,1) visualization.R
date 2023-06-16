################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)
library(xtable)
library(Hmisc)
library(kableExtra)
library(plot.matrix)
library(arfima)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

PROCESS = "FARIMA"
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14
N_SIMULATIONS = 1000

################################################################################
# LOAD data
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/",PROCESS,"/")

ar_coef_vec = readRDS(file = paste0(path,"ar_coef_vec.RData"))
ma_coef_vec = readRDS(file = paste0(path,"ma_coef_vec.RData"))
d_coef_vec = readRDS(file = paste0(path,"d_coef_vec.RData"))

################################################################################
# Running time
################################################################################

time_own = readRDS(file = paste0(path,"time_own.RData"))
time_r = readRDS(file = paste0(path,"time_R.RData"))

# Limits
lim_inf = min(c(time_own,time_r))
lim_sup = max(c(time_own,time_r))

# Two boxplots side by side
par(mfrow=c(1,2))
boxplot(time_own, ylim=c(0,0.09), names = names, ylab = "time (s)", xlab = "T")
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
boxplot(time_r, ylim=c(0,0.09), names = names, ylab = "time (s)", xlab = "T")
title(main = "farima", cex.main = 0.8, line = 0.5)
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)

graph_name = "Figure 1.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

# One graph comparing both methods
par(mfrow=c(1,1))
plot(x = POWER, y = colMeans(time_own), col = "blue", type = "o", ylim=c(lim_inf,lim_sup/2), ylab = "time (s)", labels = FALSE, xlab = "T")
lines(x = POWER, y = colMeans(time_r), col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = c("fepxmst", "farima"), col = c("blue", "red"), lty = 1)
main = paste0("Average running time for $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)
graph_name = "Figure 2.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

fit_own_phi <- readRDS(file = paste0(path, "phi_1_own.RData"))
fit_r_phi <- readRDS(file = paste0(path, "phi_1_r.RData"))

mean_phi_mse_own <- colMeans((fit_own_phi - ar_coef_vec)^2)
mean_phi_mse_r <- colMeans((fit_r_phi - ar_coef_vec)^2)

par(mfrow = c(1, 1))

# Calculate standard deviation for each column
sd_phi_mse_own <- apply(((fit_own_phi - ar_coef_vec)^2), 2, sd)
sd_phi_mse_r <- apply(((fit_r_phi - ar_coef_vec)^2), 2, sd)

# Plot mean with error bars for 'fit_own_phi'
min_lim = min(mean_phi_mse_own-sd_phi_mse_own,mean_phi_mse_r-sd_phi_mse_r)
max_lim = max(mean_phi_mse_own+sd_phi_mse_own,mean_phi_mse_r+sd_phi_mse_r)
par(mfrow = c(1, 1), mar=c(5,5,3,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_phi_mse_own, col = "blue", type = "o", ylab = TeX("$\\hat{phi}_1$ MSE"), labels = FALSE, xlab = "T", ylim = c(min_lim, max_lim))
lines(x = POWER, y = mean_phi_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col = "black")
legend("topright", legend = c("fepxmst", "farima"), col = c("blue", "red"), lty = 1)

for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i], y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] - sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] + sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")

  segments(x0 = POWER[i], y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i], y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] - sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] + sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
}
main = paste0("MSE of $\\hat{phi}_1$ with $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)
graph_name = "Figure 3.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# MSE COEFFICIENTS for theta
################################################################################

fit_own_theta = readRDS(file = paste0(path,"theta_1_own.RData"))
fit_r_theta = readRDS(file = paste0(path,"theta_1_r.RData"))

mean_theta_mse_own <- colMeans((fit_own_theta - ma_coef_vec)^2)
mean_theta_mse_r <- colMeans((fit_r_theta - ma_coef_vec)^2)

# Calculate standard deviation for each column
sd_theta_mse_own <- apply(((fit_own_theta - ma_coef_vec)^2), 2, sd)
sd_theta_mse_r <- apply(((fit_r_theta - ma_coef_vec)^2), 2, sd)

# Plot mean with error bars for 'fit_own_theta'
min_lim = min(mean_theta_mse_own-sd_theta_mse_own,mean_theta_mse_r-sd_theta_mse_r)
max_lim = max(mean_theta_mse_own+sd_theta_mse_own,mean_theta_mse_r+sd_theta_mse_r)

par(mfrow = c(1, 1), mar=c(5,5,3,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_theta_mse_own, col = "blue", type = "o", ylab = TeX("$\\hat{theta}_1$ $MSE"), labels = FALSE, xlab = "T",ylim = c(min_lim,max_lim))
lines(x = POWER, y = mean_theta_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col="black")
legend("topright", legend = c("fepxmst", "farima"), col = c("blue", "red"), lty = 1)

# Add error bars for 'fit_own_theta'
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_theta_mse_own[i] - sd_theta_mse_own[i], x1 = POWER[i], y1 = mean_theta_mse_own[i] + sd_theta_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_own[i] - sd_theta_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_own[i] - sd_theta_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_own[i] + sd_theta_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_own[i] + sd_theta_mse_own[i], col = "blue")

  segments(x0 = POWER[i], y0 = mean_theta_mse_r[i] - sd_theta_mse_r[i], x1 = POWER[i], y1 = mean_theta_mse_r[i] + sd_theta_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_r[i] - sd_theta_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_r[i] - sd_theta_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_r[i] + sd_theta_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_r[i] + sd_theta_mse_r[i], col = "red")
}
main = paste0("MSE of $\\hat{theta}_1$ with $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)
graph_name = "Figure 4.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# MSE COEFFICIENTS for d
################################################################################

fit_own_d = readRDS(file = paste0(path,"d_own.RData"))
fit_r_d = readRDS(file = paste0(path,"d_r.RData"))

mean_d_mse_own <- colMeans((fit_own_d - d_coef_vec)^2)
mean_d_mse_r <- colMeans((fit_r_d - d_coef_vec)^2)

# Calculate standard deviation for each column
sd_d_mse_own <- apply(((fit_own_d - d_coef_vec)^2), 2, sd)
sd_d_mse_r <- apply(((fit_r_d - d_coef_vec)^2), 2, sd)

# Plot mean with error bars for 'fit_own_theta'
min_lim = min(mean_d_mse_own-sd_d_mse_own,mean_d_mse_r-sd_d_mse_r)
max_lim = max(mean_d_mse_own+sd_d_mse_own,mean_d_mse_r+sd_d_mse_r)
par(mfrow = c(1, 1), mar=c(5,5,3,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_d_mse_own, col = "blue", type = "o", ylab = TeX("$\\hat{theta}_1$ $MSE"), labels = FALSE, xlab = "T",ylim = c(min_lim, max_lim))
lines(x = POWER, y = mean_d_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col="black")
legend("topright", legend = c("fepxmst", "farima"), col = c("blue", "red"), lty = 1)

# Add error bars for 'fit_own_theta'
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_d_mse_own[i] - sd_d_mse_own[i], x1 = POWER[i], y1 = mean_d_mse_own[i] + sd_d_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own[i] - sd_d_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own[i] - sd_d_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own[i] + sd_d_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own[i] + sd_d_mse_own[i], col = "blue")

  segments(x0 = POWER[i], y0 = mean_d_mse_r[i] - sd_d_mse_r[i], x1 = POWER[i], y1 = mean_d_mse_r[i] + sd_d_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_r[i] - sd_d_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_r[i] - sd_d_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_r[i] + sd_d_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_r[i] + sd_d_mse_r[i], col = "red")
}
main = paste0("MSE of $\\hat{d}$ with $\\",N_SIMULATIONS," \\ \\{y_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)
graph_name = "Figure 4.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

fit_own_exp = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp = readRDS(file = paste0(path,"lambda_r.RData"))

lim_inf = min(c(fit_own_exp,fit_r_exp))
lim_sup = max(c(fit_own_exp,fit_r_exp))

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
boxplot(fit_own_exp, ylab = TeX("$\\hat{lambda}$"), xlab = "T", ylim=c(lim_inf,lim_sup), names = names)
abline(h=1, col = "red")
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
boxplot(fit_r_exp, ylab = TeX("$\\hat{lambda}$"), xlab = "T", ylim=c(lim_inf,lim_sup), names = names)
abline(h=1, col = "red")
title(main = "farima", cex.main = 0.8, line = 0.5)
main = paste0("Fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])}(\\omega_k)_t\\}_{k=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)

graph_name = "Figure 5.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

mean_own_exp <- colMeans(fit_own_exp)
mean_r_exp <- colMeans(fit_r_exp)
sd_own_exp <- apply(fit_own_exp, 2, sd)
sd_r_exp <- apply(fit_r_exp, 2, sd)

min_lim = min(c((mean_own_exp-sd_own_exp),(mean_r_exp-sd_r_exp)))
max_lim = max(c((mean_own_exp+sd_own_exp),(mean_r_exp+sd_r_exp)))

# Plot means with error bars
par(mfrow=c(1,1), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_own_exp, col = "blue", type = "o", ylab = TeX("$\\hat{lambda}$"), xlab = "T", ylim = c(min_lim, max_lim), labels = FALSE)
lines(x = POWER, y = mean_r_exp, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)

legend("topright", legend = c("fepxmst", "farima"), col = c("blue", "red"), lty = 1)
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

main = paste0("Fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])}(\\omega_k)_t\\}_{k=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)
graph_name = "Figure 6.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

################################################################################
# Correct p-values for EXP(1)
################################################################################

p_val_own_exp = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp = readRDS(file = paste0(path,"p.val_r.RData"))

hist(p_val_own_exp)

boxplot(p_val_own_exp)
boxplot(p_val_r_exp)

lim_inf = min(c(p_val_own_exp,p_val_r_exp))
lim_sup = max(c(p_val_own_exp,p_val_r_exp))

par(mfrow=c(1,2), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
boxplot(p_val_own_exp, ylab = TeX("$p.value$"), xlab = "T", ylim=c(lim_inf,lim_sup), names = names)
abline(h=0.05, col = "red")
title(main = "fexpmst", cex.main = 0.8, line = 0.5)
boxplot(p_val_r_exp, ylab = TeX("$p.value$"), xlab = "T", ylim=c(lim_inf,lim_sup), names = names)
abline(h=0.05, col = "red")
title(main = "farima", cex.main = 0.8, line = 0.5)
main = paste0("$H_0: \\ \\{I^*_{",PROCESS,"(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\]),\\ d \\ \\sim \\ u(\\[0,0.5\\])}(\\omega_k)_t\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.2, font = 2)

graph_name = "Figure 7.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()




