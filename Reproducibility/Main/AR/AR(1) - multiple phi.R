################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)
library(xtable)
library(Hmisc)
library(kableExtra)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# SIMULATION
################################################################################

fit_own_coef_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_r_coef_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_own_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_r_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
p_val_own_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
p_val_r_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
time_own_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
time_r_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

for (param_number in 1:length(ar_coef_vec)) {
  ar_coef = ar_coef_vec[param_number]
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
      start_own = Sys.time()
      y_own = sim.farima(ar = ar_coef, T = T)
      end_own = Sys.time()
      time_own[sim,j] = end_own-start_own
      # AR EXISTING PACKAGE
      start_r = Sys.time()
      y_r = arima.sim(model = list(ar = ar_coef), n = T)
      end_r = Sys.time()
      time_r[sim,j] = end_r-start_r

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
  fit_own_coef_list[[param_number]] = fit_own_coef
  fit_r_coef_list[[param_number]] = fit_r_coef
  fit_own_exp_list[[param_number]] = fit_own_exp
  fit_r_exp_list[[param_number]] = fit_r_exp
  p_val_own_exp_list[[param_number]] = p_val_own_exp
  p_val_r_exp_list[[param_number]] = p_val_r_exp
  time_own_list[[param_number]] = time_own
  time_r_list[[param_number]] = time_r
}

################################################################################
# Running time
################################################################################


################################################################################
# MSE COEFFICIENTS for phi
################################################################################

MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_coef_list[[i]]-ar_coef_vec[i])^2)
  MSE_r[i,] = colMeans((fit_r_coef_list[[i]]-ar_coef_vec[i])^2)
}

rownames(MSE_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(MSE_own<=MSE_r)

MSE_own = xtable(MSE_own, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)

print(MSE_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_exp_list[[i]]-1)^2)
  MSE_r[i,] = colMeans((fit_r_exp_list[[i]]-1)^2)
}

rownames(MSE_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(MSE_own<=MSE_r)

MSE_own = xtable(MSE_own, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)

print(MSE_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# OWN CODE
prop_non_rejection_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
prop_non_rejection_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  prop_non_rejection_own[i,] = colMeans(p_val_own_exp_list[[i]]>=0.05)
  prop_non_rejection_r[i,] = colMeans(p_val_r_exp_list[[i]]>=0.05)
}

rownames(prop_non_rejection_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(prop_non_rejection_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(prop_non_rejection_own<prop_non_rejection_r)

prop_non_rejection_own = xtable(prop_non_rejection_own, digits = 5)
prop_non_rejection_r = xtable(prop_non_rejection_r, digits = 5)

print(prop_non_rejection_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(prop_non_rejection_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})
