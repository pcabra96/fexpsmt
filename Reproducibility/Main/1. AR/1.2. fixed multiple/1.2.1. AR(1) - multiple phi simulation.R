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

################################################################################
##----------------------------------------------------------------------------##
## 2. SEED                                                                    ##
##----------------------------------------------------------------------------##
################################################################################

set.seed(0)

################################################################################
##----------------------------------------------------------------------------##
## 3. SIMULATION PARAMETERS                                                   ##
##----------------------------------------------------------------------------##
################################################################################

PROCESS = "AR"
SUBPROCESS = "fixed multiple"
ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/1. ",PROCESS,"/1.2. ",SUBPROCESS,"/")

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION                                                              ##
##----------------------------------------------------------------------------##
################################################################################

fit_own_coef_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_r_coef_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_own_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
fit_r_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
p_val_own_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
p_val_r_exp_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
time_own_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
time_r_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
average_own_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)
average_r_list = replicate(10, matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER)), simplify = FALSE)

time_start = Sys.time()

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
  average_own = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
  average_r = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))


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
      # AVERAGE
      ################################################################################

      average_own[sim,j] = mean(y_own)
      average_r[sim,j] = mean(y_r)

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
  average_own_list[[param_number]] = average_own
  average_r_list[[param_number]] = average_r

}

time_end = Sys.time()
time = time_end-time_start
print(time)

################################################################################
##----------------------------------------------------------------------------##
## 5. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. TIME
################################################################################

saveRDS(time_own_list, file = paste0(path,"time_own_list.RData"))
saveRDS(time_r_list, file = paste0(path,"time_r_list.RData"))

################################################################################
# 5.2. AVERAGE
################################################################################

saveRDS(average_own_list, file = paste0(path,"average_own_list.RData"))
saveRDS(average_r_list, file = paste0(path,"average_r_list.RData"))

################################################################################
# 5.3. TIME DOMAIN PARAMETER AR
################################################################################

saveRDS(fit_own_coef_list, file = paste0(path,"fit_own_coef_list.RData"))
saveRDS(fit_r_coef_list, file = paste0(path,"fit_r_coef_list.RData"))

################################################################################
# 5.4. FREQUENCY DOMAIN PARAMETER
################################################################################

saveRDS(fit_own_exp_list, file = paste0(path,"fit_own_exp_list.RData"))
saveRDS(fit_r_exp_list, file = paste0(path,"fit_r_exp_list.RData"))

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

saveRDS(p_val_own_exp_list, file = paste0(path,"p_val_own_exp_list.RData"))
saveRDS(p_val_r_exp_list, file = paste0(path,"p_val_r_exp_list.RData"))

