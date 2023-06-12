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

PROCESS = "ARMA"
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14

################################################################################
# LOAD data
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/",PROCESS,"/")

################################################################################
# Running time
################################################################################

time_own = readRDS(file = paste0(path,"time_own.RData"))
time_r = readRDS(file = paste0(path,"time_R.RData"))

par(mfrow=c(1,2))
boxplot(time_own, ylim=c(0,0.09))
boxplot(time_r, ylim=c(0,0.09))

par(mfrow=c(1,1))
plot(x = POWER, y = colMeans(time_own), col = "blue", type = "l", ylim=c(0,0.03), ylab = "time (s)", labels = FALSE, xlab = "T")
lines(x = POWER, y = colMeans(time_r), col = "red", type = "l")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

fit_own_phi = readRDS(file = paste0(path,"phi_1_own.RData"))
fit_r_phi = readRDS(file = paste0(path,"phi_1_r.RData"))

fit_own_phi[,1]-

par(mfrow=c(1,1))
plot(x = POWER, y = colMeans(fit_own_phi), col = "blue", type = "l", ylim=c(0,0.03), ylab = "time (s)", labels = FALSE, xlab = "T")
lines(x = POWER, y = colMeans(fit_r_phi), col = "red", type = "l")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = c("fepxmst", "stats"), col = c("blue", "red"), lty = 1)


################################################################################
# MSE COEFFICIENTS for theta
################################################################################

fit_own_theta = readRDS(file = paste0(path,"theta_1_own.RData"))
fit_r_theta = readRDS(file = paste0(path,"theta_1_r.RData"))

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
fit_own_exp = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp = readRDS(file = paste0(path,"lambda_r.RData"))

################################################################################
# Correct p-values for EXP(1)
################################################################################

p_val_own_exp = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp = readRDS(file = paste0(path,"p.val_r.RData"))
