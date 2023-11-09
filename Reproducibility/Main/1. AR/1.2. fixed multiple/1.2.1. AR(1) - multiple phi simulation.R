################################################################################
##----------------------------------------------------------------------------##
## INDEX                                                                      ##
##----------------------------------------------------------------------------##
################################################################################

# 0. PACKAGES
# 1. SEED
# 2. SIMULATION PARAMETERS
# 3. SIMULATION
# 4. RESULTS
# 4.1. TIME
# 4.2. TIME SERIES AVERAGE
# 4.3. PHI COEFFICIENT
# 4.4. FREQUENCY DOMAIN PARAMETER
# 4.5. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
##----------------------------------------------------------------------------##
## 0. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(devtools)
devtools::install_github("pcabra96/fexpsmt", force = TRUE)
library(fexpsmt)
library(fitdistrplus)

################################################################################
##----------------------------------------------------------------------------##
## 1. SEED                                                                    ##
##----------------------------------------------------------------------------##
################################################################################

set.seed(0)

################################################################################
##----------------------------------------------------------------------------##
## 2. SIMULATION PARAMETERS                                                   ##
##----------------------------------------------------------------------------##
################################################################################

PROCESS = "AR"
SUBPROCESS = "fixed multiple"
ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
active_path = dirname(rstudioapi::getActiveDocumentContext()$path)
path = paste0(active_path,"/")

################################################################################
##----------------------------------------------------------------------------##
## 3. SIMULATION (expected runtime in mac with m1 chip: 32.74463 mins)        ##
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
## 4. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 4.1. TIME
################################################################################

saveRDS(time_own_list, file = paste0(path,"time_own_list.RData"))
saveRDS(time_r_list, file = paste0(path,"time_r_list.RData"))

################################################################################
# 4.2. TIME SERIES AVERAGE
################################################################################

saveRDS(average_own_list, file = paste0(path,"average_own_list.RData"))
saveRDS(average_r_list, file = paste0(path,"average_r_list.RData"))

################################################################################
# 4.3. PHI COEFFICIENT
################################################################################

saveRDS(fit_own_coef_list, file = paste0(path,"fit_own_coef_list.RData"))
saveRDS(fit_r_coef_list, file = paste0(path,"fit_r_coef_list.RData"))

################################################################################
# 4.4. FREQUENCY DOMAIN PARAMETER
################################################################################

saveRDS(fit_own_exp_list, file = paste0(path,"fit_own_exp_list.RData"))
saveRDS(fit_r_exp_list, file = paste0(path,"fit_r_exp_list.RData"))

################################################################################
# 4.5. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

saveRDS(p_val_own_exp_list, file = paste0(path,"p_val_own_exp_list.RData"))
saveRDS(p_val_r_exp_list, file = paste0(path,"p_val_r_exp_list.RData"))

