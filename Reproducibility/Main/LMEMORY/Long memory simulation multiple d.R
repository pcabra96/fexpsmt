################################################################################
# INDEX
################################################################################

# 1. PACKAGES
# 2. SEED
# 3. SIMULATION PARAMETERS
# 4. SIMULATION
# 5. RESULTS
# 5.1. TIME
# 5.2. TIME DOMAIN PARAMETER
# 5.3. FREQUENCY DOMAIN PARAMETER
# 5.3. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
# 1. PACKAGES
################################################################################

library(arfima)
library(forecast)
require(MASS)
library(latex2exp)

################################################################################
# 2. SEED
################################################################################

set.seed(0)

################################################################################
# 3. SIMULATION PARAMETERS
################################################################################

PROCESS = "FARIMA"
symbol = "d"
d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

# TIME
own_times_farima_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
own_times_fexp_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
r_times_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

# PARAMETER ESTIMATION TIME DOMAIN
own_long_param_farima_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
own_long_param_fexp_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
r_long_param_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

# PARAMETER ESTIMATION FREQUENCY DOMAIN
own_lambda_farima_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
own_lambda_fexp_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
r_lambda_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

# GOODNESS OF FIT
own_exp_1_farima_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
own_exp_1_fexp_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
r_exp_1_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

################################################################################
# 4. SIMULATION
################################################################################
begin_time = Sys.time()
for (param_number in 1:length(d_coef_vec)) {

  d_coef = d_coef_vec[param_number]

  # TIME
  own_times_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  own_times_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  r_times = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

  # PARAMETER ESTIMATION TIME DOMAIN
  own_long_param_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  own_long_param_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  r_long_param = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

  # PARAMETER ESTIMATION FREQUENCY DOMAIN
  own_lambda_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  own_lambda_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  r_lambda = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

  # GOODNESS OF FIT
  own_exp_1_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  own_exp_1_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  r_exp_1 = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

  for (sim in 1:N_SIMULATIONS) {
    for (j in 1:length(POWER)) {
      T = 2^POWER[j]
      true_spectrum = farima.spectrum(d = d_coef, n.freq = T)
      true_spectrum[1] = true_spectrum[length(true_spectrum)]

      # OWN TIME FOR FARIMA
      own_start = Sys.time()
      y_own_farima  = sim.farima(d = d_coef, T = T)
      own_end = Sys.time()
      own_times_farima[sim,j] = own_end-own_start

      # OWN TIME FOR FARIMA
      own_start = Sys.time()
      y_own_fexp  = sim.fexp(ck = 0, d = d_coef, T = T)
      own_end = Sys.time()
      own_times_fexp[sim,j] = own_end-own_start

      # R TIME
      r_start = Sys.time()
      y_R = arfima.sim(n = T, model = list(dfrac=d_coef))
      r_end = Sys.time()
      r_times[sim,j] = r_end-r_start

      ##################################
      # FITTING
      ##################################

      own_long_param_farima[sim,j] = fit.farima(y_own_farima, d = 1)[["d"]]
      own_long_param_fexp[sim,j] = fit.farima(y_own_fexp, d = 1)[["d"]]
      r_long_param[sim,j] = fit.farima(y_R, d = 1)[["d"]]

      ################################################################################
      # PERIODOGRAM
      ################################################################################

      n = length(y_own_farima)

      # Fundamental frequencies
      mhalfm <- (n-1) %/% 2L
      w <- 2*pi/n * (1:mhalfm)

      # Periodogram ordinates by FFT of own simulation for FARIMA
      per_own_farima = (Mod(fft(y_own_farima))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
      yper_own_farima = per_own_farima[2: ((n+1) %/% 2)]
      I_own_farima = yper_own_farima/true_spectrum[1:(T/2-1)]

      # Periodogram ordinates by FFT of own simulation for FEXP
      per_own_fexp = (Mod(fft(y_own_fexp))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
      yper_own_fexp = per_own_fexp[2: ((n+1) %/% 2)]
      I_own_fexp = yper_own_fexp/true_spectrum[1:(T/2-1)]

      # Periodogram ordinates by FFT of R simulation
      per_r = (Mod(fft(y_R))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
      yper_r = per_r[2: ((n+1) %/% 2)]
      I_r = yper_r/true_spectrum[1:(T/2-1)]

      # Periodogram distribution of own simulation for FARIMA
      own_lambda_farima[sim,j] <- fitdistr(I_own_farima, "exponential")[["estimate"]]
      own_exp_1_farima[sim,j] = ks.test(I_own_farima, "pexp", 1)[["p.value"]]

      # Periodogram distribution of own simulation for FEXP
      own_lambda_fexp[sim,j] <- fitdistr(I_own_fexp, "exponential")[["estimate"]]
      own_exp_1_fexp[sim,j] = ks.test(I_own_fexp, "pexp", 1)[["p.value"]]

      # Periodogram distribution of R simulation
      r_lambda[sim,j] <- fitdistr(I_r, "exponential")[["estimate"]]
      r_exp_1[sim,j] = ks.test(I_r, "pexp", 1)[["p.value"]]
    }
  }

  # TIME
  own_times_farima_list[[param_number]] = own_times_farima
  own_times_fexp_list[[param_number]] = own_times_fexp
  r_times_list[[param_number]] = r_times

  # PARAMETER ESTIMATION TIME DOMAIN
  own_long_param_farima_list[[param_number]] = own_long_param_farima
  own_long_param_fexp_list[[param_number]] = own_long_param_fexp
  r_long_param_list[[param_number]] = r_long_param

  # PARAMETER ESTIMATION FREQUENCY DOMAIN
  own_lambda_farima_list[[param_number]] = own_lambda_farima
  own_lambda_fexp_list[[param_number]] = own_lambda_fexp
  r_lambda_list[[param_number]] = r_lambda

  # GOODNESS OF FIT
  own_exp_1_farima_list[[param_number]] = own_exp_1_farima
  own_exp_1_fexp_list[[param_number]] = own_exp_1_fexp
  r_exp_1_list[[param_number]] = r_exp_1
}
end_time = Sys.time()
total_time = end_time - begin_time
################################################################################
# VISUALIZATION
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/LMEMORY/")

########################################################
# 5.1. TIME
########################################################

# Lines
min_lim = min(r_times,own_times_farima,own_times_fexp)
max_lim = max(r_times,own_times_farima,own_times_fexp)

# Boxplot
par(mfrow=c(1,3), mar=c(5,5,5,2)) # mar = c(bottom, left, top, right))
min_lim = min(r_times,own_times_farima,own_times_fexp)
max_lim = max(r_times,own_times_farima,own_times_fexp)
boxplot(own_times_farima, ylim=c(min_lim, max_lim), names = names, ylab = "time (s)", xlab = "T")
title(main = TeX(paste0("$fexpmst_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
boxplot(own_times_fexp, ylim=c(min_lim, max_lim), names = names, ylab = "time (s)", xlab = "T")
title(main = TeX(paste0("$fexpmst_{FEXP(0,",d_coef,")}$")), cex.main = 1, line = 0.5)
boxplot(r_times, ylim=c(min_lim, max_lim), names = names, ylab = "time (s)", xlab = "T")
title(main = TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
main = paste0("Simulation time for $\\",N_SIMULATIONS,"\\{ y_{FARIMA(0,",d_coef,",0)_t} \\ or \\ y_{FEXP(0,",d_coef,"))_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -3, outer = TRUE,cex=1.4, font = 2)

graph_name = "Figure 1.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

# Lines
min_lim = min(colMeans(r_times),colMeans(own_times_farima),colMeans(own_times_fexp))
max_lim = max(colMeans(r_times),colMeans(own_times_farima),colMeans(own_times_fexp))

par(mfrow=c(1,1))
plot(x = POWER, y = colMeans(own_times_farima), type = "o", ylim=c(min_lim,max_lim), col = "blue", ylab = "time (s)", xlab = "T", labels = FALSE)
lines(x = POWER, y = colMeans(own_times_fexp), type = "o", col = "green")
lines(x = POWER, y = colMeans(r_times), type = "o", col = "red")
legend("topleft", legend = c(TeX(paste0("$fepxmst_{FARIMA(0,",d_coef,",0)}$")),TeX(paste0("$fepxmst_{FEXP(0,",d_coef,",0)}$")), TeX(paste0("$arfima_{FARIMA(0,",d_coef,",0)}$"))), col = c("blue", "green", "red"), lty = 1)
axis(1, at=POWER, labels = names)
axis(2)
main = paste0("Average running time for $\\",N_SIMULATIONS,"\\ \\{ y_{FARIMA(0,",d_coef,",0)_t \\ or \\ FEXP(0,",d_coef,"))_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 2.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

########################################################
# 5.2. TIME DOMAIN PARAMETER PARAMETER
########################################################

# Fitted d coefficient BOXPLOTS
par(mfrow=c(1,3), mar=c(5,5,5,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{",symbol,"}$")
boxplot(own_long_param_farima, names = names, xlab = "T", ylab=TeX(main))
title(main = TeX(paste0("$fexpmst_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
abline(h=d_coef, col = "red")
main = paste0("$\\hat{",symbol,"}$")
boxplot(own_long_param_fexp, names = names, xlab = "T", ylab = TeX(main))
title(main = TeX(paste0("$fexpmst_{FEXP(0,",d_coef,")}$")), cex.main = 0.8, line = 0.5)
abline(h=d_coef, col = "red")
main = paste0("$\\hat{",symbol,"}$")
boxplot(r_long_param, names = names, xlab = "T", ylab = TeX(main))
title(main = TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$")), cex.main = 0.8, line = 0.5)
abline(h=d_coef, col = "red")

main = paste0("Fitted $",symbol,"$ for$\\ ",N_SIMULATIONS," \\ \\{ y_{FARIMA(0,",d_coef,",0)_t \\ or \\ FEXP(0,",d_coef,"))_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -3, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 3.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)

# Fitted d coefficient LINES

# Calculate standard deviation for each column
mean_d_mse_own_farima <- colMeans((own_long_param_farima - d_coef)^2)
mean_d_mse_own_fexp <- colMeans((own_long_param_fexp - d_coef)^2)
mean_d_mse_r <- colMeans((r_long_param - d_coef)^2)
sd_d_mse_own_farima <- apply(((own_long_param_farima - d_coef)^2), 2, sd)
sd_d_mse_own_fexp <- apply(((own_long_param_fexp - d_coef)^2), 2, sd)
sd_d_mse_r <- apply(((r_long_param - d_coef)^2), 2, sd)

# Plot mean with error bars for 'fit_own_phi'
min_lim = min(mean_d_mse_own_farima-sd_d_mse_own_farima,mean_d_mse_own_fexp-sd_d_mse_own_fexp,mean_d_mse_r-sd_d_mse_r)
max_lim = max(mean_d_mse_own_farima+sd_d_mse_own_farima,mean_d_mse_own_fexp+sd_d_mse_own_fexp,mean_d_mse_r+sd_d_mse_r)

par(mfrow = c(1, 1), mar=c(5,5,3,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{",symbol,"}$ MSE")
plot(x = POWER, y = mean_d_mse_own_farima, col = "blue", type = "o", ylab = TeX(main), labels = FALSE, xlab = "T", ylim = c(min_lim,max_lim))
lines(x = POWER, y = mean_d_mse_own_fexp, col = "green", type = "o")
lines(x = POWER, y = mean_d_mse_r, col = "red", type = "o")
axis(1, at = POWER, labels = names)
axis(2)
abline(h=0,col = "black")
legend("topright", legend = c(TeX(paste0("$fepxmst_{FARIMA(0,",d_coef,",0)}$")),TeX(paste0("$fepxmst_{FEXP(0,",d_coef,")}$")), TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$"))), col = c("blue", "green", "red"), lty = 1)

for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_d_mse_own_farima[i] - sd_d_mse_own_farima[i], x1 = POWER[i], y1 = mean_d_mse_own_farima[i] + sd_d_mse_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own_farima[i] - sd_d_mse_own_farima[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own_farima[i] - sd_d_mse_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own_farima[i] + sd_d_mse_own_farima[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own_farima[i] + sd_d_mse_own_farima[i], col = "blue")

  segments(x0 = POWER[i], y0 = mean_d_mse_own_fexp[i] - sd_d_mse_own_fexp[i], x1 = POWER[i], y1 = mean_d_mse_own_fexp[i] + sd_d_mse_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own_fexp[i] - sd_d_mse_own_fexp[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own_fexp[i] - sd_d_mse_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_own_fexp[i] + sd_d_mse_own_fexp[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_own_fexp[i] + sd_d_mse_own_fexp[i], col = "green")

  segments(x0 = POWER[i], y0 = mean_d_mse_r[i] - sd_d_mse_r[i], x1 = POWER[i], y1 = mean_d_mse_r[i] + sd_d_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_r[i] - sd_d_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_r[i] - sd_d_mse_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_d_mse_r[i] + sd_d_mse_r[i], x1 = POWER[i] + 0.1, y1 = mean_d_mse_r[i] + sd_d_mse_r[i], col = "red")
}
main = paste0("MSE of $\\hat{",symbol,"}$ with $\\",N_SIMULATIONS,"\\ \\{y_{FARIMA(0,",d_coef,",0)_t} \\ or \\ y_{FEXP(0,",d_coef,"))_t}\\}_{t=1}^{T}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 4.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)

########################################################
# 5.3. FREQUENCY DOMAIN PARAMETER PARAMETER
########################################################

# Fitted FARIMA and FEXP periodogram BOXPLOTS
min_lim = min(own_lambda_farima, own_lambda_fexp, r_lambda)
max_lim = max(own_lambda_farima, own_lambda_fexp, r_lambda)

par(mfrow = c(1, 3), mar=c(5,5,5,2)) # mar = c(bottom, left, top, right))
main = paste0("$\\hat{lambda}_{MLE}$")
boxplot(own_lambda_farima, names = names, xlab = "T", ylab = TeX(main), ylim=c(min_lim,max_lim))
title(main = TeX(paste0("$fexpmst_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
abline(h=1, col = "red")
main = paste0("$\\hat{lambda}_{MLE}$")
boxplot(own_lambda_fexp, names = names, xlab = "T", ylab =TeX(main), ylim=c(min_lim,max_lim))
title(main = TeX(paste0("$fexpmst_{FEXP(0,",d_coef,")}$")), cex.main = 1, line = 0.5)
abline(h=1, col = "red")
main = paste0("$\\hat{lambda}_{MLE}$")
boxplot(r_lambda, names = names, xlab = "T", ylab =TeX(main), ylim=c(min_lim,max_lim))
title(main = TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
abline(h=1, col = "red")
main = paste0("Fitted $\\lambda$ for ", N_SIMULATIONS, " $\\{I^*_{FARIMA(0,",d_coef,",0)_t} \\ or \\ I^*_{FEXP(0,",d_coef,"))_t}\\}_{t=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 5.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)

# Fitted FARIMA and FEXP periodogram LINES

mean_own_exp_farima <- colMeans((own_lambda_farima-1)^2)
mean_own_exp_fexp <- colMeans((own_lambda_fexp-1)^2)
mean_r_exp <- colMeans((r_lambda-1)^2)

sd_own_exp_farima <- apply((own_lambda_farima-1)^2, 2, sd)
sd_own_exp_fexp <- apply((own_lambda_fexp-1)^2, 2, sd)
sd_r_exp <- apply((r_lambda-1)^2, 2, sd)

min_lim = min(c((mean_own_exp_farima-sd_own_exp_farima),(mean_own_exp_fexp-sd_own_exp_fexp),(mean_r_exp-sd_r_exp)))
max_lim = max(c((mean_own_exp_farima+sd_own_exp_farima),(mean_own_exp_fexp+sd_own_exp_fexp),(mean_r_exp+sd_r_exp)))

par(mfrow=c(1,1), mar=c(5,5,4,2)) # mar = c(bottom, left, top, right))
plot(x = POWER, y = mean_own_exp_farima, col = "blue", type = "o", ylab = TeX("MSE $\\hat{lambda}$"), xlab = "T", ylim = c(min_lim, max_lim), labels = FALSE)
lines(x = POWER, y = mean_own_exp_fexp, col = "green", type = "o")
lines(x = POWER, y = mean_r_exp, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)

legend("topright", legend = c(TeX(paste0("$fepxmst_{FARIMA(0,",d_coef,",0)}$")), TeX(paste0("$fepxmst_{FEXP(0,",d_coef,")}$")),TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$"))), col = c("blue", "green","red"), lty = 1)
abline(h=1, col = "black")

# Add error bars
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = mean_own_exp_farima[i] - sd_own_exp_farima[i], x1 = POWER[i], y1 = mean_own_exp_farima[i] + sd_own_exp_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp_farima[i] - sd_own_exp_farima[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp_farima[i] - sd_own_exp_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp_farima[i] + sd_own_exp_farima[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp_farima[i] + sd_own_exp_farima[i], col = "blue")

  segments(x0 = POWER[i], y0 = mean_own_exp_fexp[i] - sd_own_exp_fexp[i], x1 = POWER[i], y1 = mean_own_exp_fexp[i] + sd_own_exp_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp_fexp[i] - sd_own_exp_fexp[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp_fexp[i] - sd_own_exp_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = mean_own_exp_fexp[i] + sd_own_exp_fexp[i], x1 = POWER[i] + 0.1, y1 = mean_own_exp_fexp[i] + sd_own_exp_fexp[i], col = "green")

  segments(x0 = POWER[i], y0 = mean_r_exp[i] - sd_r_exp[i], x1 = POWER[i], y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_r_exp[i] - sd_r_exp[i], x1 = POWER[i] + 0.1, y1 = mean_r_exp[i] - sd_r_exp[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = mean_r_exp[i] + sd_r_exp[i], x1 = POWER[i] + 0.1, y1 = mean_r_exp[i] + sd_r_exp[i], col = "red")
}

main = paste0("MSE of fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{FARIMA(0,",d_coef,",0)_t} or  I^*_{FEXP(0,",d_coef,")_t}\\}_{k=1}^{T-1}$")
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 6.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()

########################################################
# 5.4. GOODNESS OF FIT FREQUENCY DOMAIN
########################################################

# P value
par(mfrow = c(1, 3), mar=c(5,5,5,2)) # mar = c(bottom, left, top, right))
boxplot(own_exp_1_farima, ylab = "p.value", names = names, xlab = "T")
title(main = TeX(paste0("$fexpmst_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
abline(h=0.05, col = "red")

boxplot(own_exp_1_fexp, ylab = "p.value", names = names, xlab = "T")
title(main = TeX(paste0("$fexpmst_{FEXP(0,",d_coef,")}$")), cex.main = 1, line = 0.5)
abline(h=0.05, col = "red")

boxplot(r_exp_1, ylab = "p.value", names = names, xlab = "T")
title(main = TeX(paste0("$farima_{FARIMA(0,",d_coef,",0)}$")), cex.main = 1, line = 0.5)
abline(h=0.05, col = "red")

main = paste0("$H_0: \\ \\{I^*_{FARIMA(0,",d_coef,",0)_t} \\ or \\ I^*_{FEXP(0,",d_coef,")_t}\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 7.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()




