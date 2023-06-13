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

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

PROCESS = "ARMA"
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14
N_SIMULATIONS = 1000

################################################################################
# LOAD data
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/",PROCESS,"/")

ar_coef_vec = readRDS(file = paste0(path,"ar_coef_vec.RData"))
ma_coef_vec = readRDS(file = paste0(path,"ma_coef_vec.RData"))

################################################################################
# Running time
################################################################################

time_own = readRDS(file = paste0(path,"time_own.RData"))
time_r = readRDS(file = paste0(path,"time_R.RData"))

par(mfrow=c(1,2))

boxplot(time_own, ylim=c(0,0.09), names = names, ylab = "time (s)", xlab = "T", main = "fexpmst")
boxplot(time_r, ylim=c(0,0.09), names = names, ylab = "time (s)", xlab = "T", main = "stats")
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{ARMA(phi_1,theta_1 \\ \\sim \\ u(\\[-0.9,0.9\\])_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -1.8, outer = TRUE,cex=1.2, )

par(mfrow=c(1,1))
plot(x = POWER, y = colMeans(time_own), col = "blue", type = "o", ylim=c(0,0.01), ylab = "time (s)", labels = FALSE, xlab = "T")
lines(x = POWER, y = colMeans(time_r), col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)
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
plot(x = POWER, y = mean_phi_mse_own, col = "blue", type = "o", ylab = TeX("$\\hat{phi}_1$ MSE"), labels = FALSE, xlab = "T", ylim = c(min(mean_phi_mse_own-sd_phi_mse_own),max(mean_phi_mse_own+sd_phi_mse_own)))
lines(x = POWER, y = mean_phi_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col = "black")
legend("topright", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)

for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i], y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] - sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] - sd_phi_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_own[i] + sd_phi_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_own[i] + sd_phi_mse_own[i], col = "blue")
}

# Add error bars for 'fit_r_phi'
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i], y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] - sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] - sd_phi_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_phi_mse_r[i] + sd_phi_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_phi_mse_r[i] + sd_phi_mse_r[i], col = "red")
}

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
plot(x = POWER, y = mean_theta_mse_own, col = "blue", type = "l", ylab = TeX("$\\hat{theta}_1$ $MSE"), labels = FALSE, xlab = "T",ylim = c(min(mean_theta_mse_own-sd_theta_mse_own),max(mean_theta_mse_own+sd_theta_mse_own)))
lines(x = POWER, y = mean_theta_mse_r, col = "red", type = "l")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col="black")
legend("topright", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)

# Add error bars for 'fit_own_theta'
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_theta_mse_own[i] - sd_theta_mse_own[i], x1 = POWER[i], y1 = mean_theta_mse_own[i] + sd_theta_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_own[i] - sd_theta_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_own[i] - sd_theta_mse_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_own[i] + sd_theta_mse_own[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_own[i] + sd_theta_mse_own[i], col = "blue")
}

# Add error bars for 'fit_r_theta'
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_theta_mse_r[i] - sd_theta_mse_r[i], x1 = POWER[i], y1 = mean_theta_mse_r[i] + sd_theta_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_r[i] - sd_theta_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_r[i] - sd_theta_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_theta_mse_r[i] + sd_theta_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_theta_mse_r[i] + sd_theta_mse_r[i], col = "red")
}
dev.off()

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

main = "$\\thta_1 \sim u(-0.9,0.9)$"

# OWN CODE
fit_own_exp = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp = readRDS(file = paste0(path,"lambda_r.RData"))

par(mfrow=c(1,2))
boxplot(fit_own_exp)
boxplot(fit_r_exp)

mean_own_exp <- colMeans(fit_own_exp)
mean_r_exp <- colMeans(fit_r_exp)
sd_own_exp <- apply(fit_own_exp, 2, sd)
sd_r_exp <- apply(fit_r_exp, 2, sd)

# Set up the x-axis labels
x_labels <- 1:length(mean_own_exp)

min_lim = min(c((mean_own_exp-sd_own_exp),(mean_r_exp-sd_r_exp)))
max_lim = max(c((mean_own_exp+sd_own_exp),(mean_r_exp+sd_r_exp)))

# Plot means with error bars
plot(x = x_labels, y = mean_own_exp, col = "blue", type = "l", ylab = "MSE", xlab = "Column", ylim = c(min_lim, max_lim))
lines(x = x_labels, y = mean_r_exp, col = "red", type = "l")
axis(1, at = x_labels, labels = x_labels)
axis(2)
legend("topright", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)
abline(h=1, col = "black")

# Add error bars
for (i in x_labels) {
  segments(x0 = i, y0 = mean_own_exp[i] - sd_own_exp[i], x1 = i, y1 = mean_own_exp[i] + sd_own_exp[i], col = "blue")
  segments(x0 = i - 0.1, y0 = mean_own_exp[i] - sd_own_exp[i], x1 = i + 0.1, y1 = mean_own_exp[i] - sd_own_exp[i], col = "blue")
  segments(x0 = i - 0.1, y0 = mean_own_exp[i] + sd_own_exp[i], x1 = i + 0.1, y1 = mean_own_exp[i] + sd_own_exp[i], col = "blue")

  segments(x0 = i, y0 = mean_r_exp[i] - sd_r_exp[i], x1 = i, y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
  segments(x0 = i - 0.1, y0 = mean_r_exp[i] - sd_r_exp[i], x1 = i + 0.1, y1 = mean_r_exp[i] - sd_r_exp[i], col = "red")
  segments(x0 = i - 0.1, y0 = mean_r_exp[i] + sd_r_exp[i], x1 = i + 0.1, y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
}

################################################################################
# Correct p-values for EXP(1)
################################################################################

p_val_own_exp = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp = readRDS(file = paste0(path,"p.val_r.RData"))

boxplot(p_val_own_exp)
boxplot(p_val_r_exp)
