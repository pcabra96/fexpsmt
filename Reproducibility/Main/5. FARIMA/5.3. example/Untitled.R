################################################################################
# PACKAGES
################################################################################

library(fracdiff)
library(latex2exp)
library(MASS)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

PROCESS = "FARIMA"
SUBPROCESS = "example"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/5. ",PROCESS,"/5.3. ",SUBPROCESS,"/")

active_path = dirname(rstudioapi::getActiveDocumentContext()$path)
path = paste0(active_path,"/")

POWER = 13
ar_coef_vec = 0.7
ma_coef_vec = 0.4
d_coef_vec = 0.25

################################################################################
# Save data
################################################################################

# TIME
time_own = 0
time_r = 0

# AVERAGE
average_own = 0
average_r = 0

# AR COEF
fit_own_phi = 0
fit_own_theta = 0

# MA COEF
fit_r_phi = 0
fit_r_theta = 0

# d coefficient
fit_r_d = 0
fit_own_d = 0

# LAMBDA
fit_own_exp = 0
fit_r_exp = 0

# P-VAL
p_val_own_exp = 0
p_val_r_exp = 0

################################################################################
# SIMULATION
################################################################################

T = 2^POWER

ar_coef = ar_coef_vec
ma_coef = ma_coef_vec
d_coef = d_coef_vec
true_spectrum = farima.spectrum(ar = ar_coef,ma= ma_coef, d = d_coef, n.freq = T)
true_spectrum[1] = true_spectrum[length(true_spectrum)]

# ARMA(1,d,1) OWN simulations
start_own = Sys.time()
y_own = sim.farima(ar = ar_coef, ma = ma_coef, d = d_coef, T = T)
end_own = Sys.time()
time_own = end_own-start_own

# AR EXISTING PACKAGE
start_r = Sys.time()
y_r = fracdiff.sim(n = T, ar = ar_coef, ma = -ma_coef, d = d_coef)[["series"]]
end_r = Sys.time()
time_r = end_r-start_r

lab_size = 1.8
axis_text = 1.5

graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow=c(2,1), mar=c(5,5,1,1))
plot(y_own, type = "l", ylab = TeX("$y_{fexpmst}$"), xlab = "T", main = "", cex.lab = lab_size, cex.axis = axis_text)
plot(y_r, type = "l", ylab = TeX("$y_{fracdiff}$"), xlab = "T", cex.lab = lab_size, cex.axis = axis_text)

dev.off()

################################################################################
# AVERAGE
################################################################################

average_own = mean(y_own)
average_r = mean(y_r)

################################################################################
# FITTING
################################################################################

a = fracdiff(x = y_own,nar = 1, nma = 1)
fit_own_phi = a[["ar"]]
fit_own_theta = -a[["ma"]]
fit_own_d = a[["d"]]

a = fracdiff(x = y_r,nar = 1, nma = 1)
fit_r_phi = a[["ar"]]
fit_r_theta = -a[["ma"]]
fit_r_d = a[["d"]]

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


graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow=c(2,1), mar=c(5,5,1,1))
plot(x = w, y = yper_own, type = "l", col = "black", ylim=c(0,max(yper_own,yper_r)), ylab = TeX("$I_{fexpmst}$"), xlab = TeX("$\\omega$"), cex.lab = lab_size, cex.axis = axis_text)
lines(x = w, y = true_spectrum[1:(T/2-1)], type = "l", col = "red")
legend("topright",legend = c("periodogram from fexpmst simulation", "true value"), col = c("black", "red"), lty = 1, cex = 1)

plot(x = w, y = yper_r, type = "l", col = "black", ylim=c(0,max(yper_own,yper_r)), ylab = TeX("$I_{fracdiff}$"), xlab = TeX("$\\omega$"), cex.lab = lab_size, cex.axis = axis_text)
lines(x = w, y = true_spectrum[1:(T/2-1)], type = "l", col = "red")
legend("topright",legend = c("periodogram from fraciff simulation", "true value"), col = c("black", "red"), lty = 1, cex = 1)

dev.off()

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow=c(2,1), mar=c(5,5,1,1))
hist(I_own, xlab = TeX("$I^*_{fexpmst}(\\omega_j)"), main = "", probability = TRUE, breaks = 25, xlim= c(0,max(I_r,I_own)), cex.lab = lab_size, cex.axis = axis_text)
curve(dexp(x, rate = 1), add = TRUE, col = "red", lwd = 2)
legend("topright",legend = c(TeX("$EXP(\\lambda=1)$")), col = "red", lty = 1, cex = 1)

hist(I_r, xlab = TeX("$I^*_{fracdiff}(\\omega_j)"), main = "", probability = TRUE, breaks = 25, xlim= c(0,max(I_r,I_own)), cex.lab = lab_size, cex.axis = axis_text)
curve(dexp(x, rate = 1), add = TRUE, col = "red", lwd = 2)
legend("topright",legend = c(TeX("$EXP(\\lambda=1)$")), col = "red", lty = 1, cex = 1)

dev.off()

# Periodogram distribution of own simulation
fit_own_exp <- fitdistr(I_own, "exponential")[["estimate"]]
p_val_own_exp = ks.test(I_own, "pexp", 1)[["p.value"]]

# Periodogram distribution of R simulation
fit_r_exp <- fitdistr(I_r, "exponential")[["estimate"]]
p_val_r_exp = ks.test(I_r, "pexp", 1)[["p.value"]]

################################################################################
# SAVE AR and MA parameteres
################################################################################

saveRDS(ar_coef_vec, file = paste0(path,"ar_coef_vec.RData"))
saveRDS(ma_coef_vec, file = paste0(path,"ma_coef_vec.RData"))
saveRDS(d_coef_vec, file = paste0(path,"d_coef_vec.RData"))

################################################################################
# Running time
################################################################################

colnames(time_own) = POWER
colnames(time_r) = POWER

saveRDS(time_own, file = paste0(path,"time_own.RData"))
saveRDS(time_r, file = paste0(path,"time_R.RData"))

################################################################################
# AVERAGE
################################################################################

colnames(average_own) = POWER
colnames(average_r) = POWER

saveRDS(average_own, file = paste0(path,"average_own.RData"))
saveRDS(average_r, file = paste0(path,"average_r.RData"))

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

colnames(fit_own_phi) = POWER
colnames(fit_r_phi) = POWER

saveRDS(fit_own_phi, file = paste0(path,"phi_1_own.RData"))
saveRDS(fit_r_phi, file = paste0(path,"phi_1_r.RData"))

################################################################################
# MSE COEFFICIENTS for theta
################################################################################

colnames(fit_own_theta) = POWER
colnames(fit_r_theta) = POWER

saveRDS(fit_own_theta, file = paste0(path,"theta_1_own.RData"))
saveRDS(fit_r_theta, file = paste0(path,"theta_1_r.RData"))

################################################################################
# MSE COEFFICIENTS for d
################################################################################

colnames(fit_own_d) = POWER
colnames(fit_r_d) = POWER

saveRDS(fit_own_d, file = paste0(path,"d_own.RData"))
saveRDS(fit_r_d, file = paste0(path,"d_r.RData"))

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
colnames(fit_own_exp) = POWER
colnames(fit_r_exp) = POWER

saveRDS(fit_own_exp, file = paste0(path,"lambda_own.RData"))
saveRDS(fit_r_exp, file = paste0(path,"lambda_r.RData"))

################################################################################
# Correct p-values for EXP(1)
################################################################################

# OWN CODE
colnames(p_val_own_exp) = POWER
colnames(p_val_r_exp) = POWER

saveRDS(p_val_own_exp, file = paste0(path,"p.val_own.RData"))
saveRDS(p_val_r_exp, file = paste0(path,"p.val_r.RData"))
