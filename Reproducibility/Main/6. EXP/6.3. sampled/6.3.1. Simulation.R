
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
# 5.2. AVERAGE
# 5.3. TIME DOMAIN PARAMETER
# 5.4. FREQUENCY DOMAIN PARAMETER
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(fracdiff)
library(latex2exp)
library(MASS)

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

PROCESS = "EXP"
SUBPROCESS = "sampled"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/6. ",PROCESS,"/6.3. ",SUBPROCESS,"/")
POWER = 7:14
N_SIM = 1000
K = 4

################################################################################
# Save data
################################################################################

# TIME
time_ar = matrix(0,nrow = N_SIM, ncol = length(POWER))
time_exp = matrix(0,nrow = N_SIM, ncol = length(POWER))

# AVERAGE
ar_average  = matrix(0,nrow = N_SIM, ncol = length(POWER))
exp_average  = matrix(0,nrow = N_SIM, ncol = length(POWER))

# PARAMETER ESTIMATION TIME DOMAIN
true_parameteres = replicate(length(POWER),matrix(0,nrow = N_SIM, ncol = K), simplify = FALSE)
true_ar = matrix(0,nrow = N_SIM, ncol = length(POWER))
fit_ar_ar = matrix(0,nrow = N_SIM, ncol = length(POWER))
fit_ar_exp = matrix(0,nrow = N_SIM, ncol = length(POWER))
fit_c_k_ar = replicate(length(POWER),matrix(0,nrow = N_SIM, ncol = K), simplify = FALSE)
fit_c_k_exp = replicate(length(POWER),matrix(0,nrow = N_SIM, ncol = K), simplify = FALSE)

# PARAMETER ESTIMATION FREQUENCY DOMAIN
ar_lambda  = matrix(0,nrow = N_SIM, ncol = length(POWER))
ar_exp_1  = matrix(0,nrow = N_SIM, ncol = length(POWER))

# GOODNESS OF FIT
exp_lambda  = matrix(0,nrow = N_SIM, ncol = length(POWER))
exp_exp_1  = matrix(0,nrow = N_SIM, ncol = length(POWER))

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION                                                              ##
##----------------------------------------------------------------------------##
################################################################################

time_start = Sys.time()
for (sim in 1:N_SIM) {
    for (j in 1:length(POWER)){
      T = 2^POWER[j]
      ar_coef = runif(n = 1, min = -0.7, max = 0.7)
      f_t = farima.spectrum(ar = ar_coef, n.freq = T)
      coef_ck = fourier.series(f_t = log(f_t), k = K)$coef
      true_ar[sim,j] = ar_coef
      true_parameteres[[j]][sim,] = coef_ck
      true_spectrum = fexp.spectrum(ck = coef_ck,n.freq = T)

      ##################################
      # FITTING
      ##################################

      t_ar_begin = Sys.time()
      y_ar = sim.farima(ar = ar_coef, T = T)
      t_ar_end = Sys.time()
      time_ar[sim,j] = abs(t_ar_end-t_ar_begin)

      t_exp_begin = Sys.time()
      y_exp = sim.fexp(ck = coef_ck, T = T)
      t_exp_end = Sys.time()
      time_exp[sim,j] = abs(t_exp_end-t_exp_begin)

      ##################################
      # FITTING
      ##################################

      fit_ar_ar[sim,j] = fit.farima(y_ar,p = 1)$ar
      fit_ar_exp[sim,j] = fit.farima(y_exp,p = 1)$ar

      fit_c_k_ar[[j]][sim,] = fit.fexp(y_ar,p = K)$c_k
      fit_c_k_exp[[j]][sim,] = fit.fexp(y_exp,p = K)$c_k

      ################################################################################
      # AVERAGE
      ################################################################################

      ar_average[sim,j] = mean(y_ar)
      exp_average[sim,j] = mean(y_exp)

      ################################################################################
      # PERIODOGRAM
      ################################################################################

      n = length(y_ar)

      # Fundamental frequencies
      mhalfm <- (n-1) %/% 2L
      w <- 2*pi/n * (1:mhalfm)

      # Periodogram ordinates by FFT of AR(1) benchmark
      per_ar = (Mod(fft(y_ar))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
      yper_ar = per_ar[2: ((n+1) %/% 2)]
      I_ar = yper_ar/true_spectrum[1:(T/2-1)]

      # Periodogram ordinates by FFT of AR(1) benchmark
      per_exp = (Mod(fft(y_exp))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
      yper_exp = per_exp[2: ((n+1) %/% 2)]
      I_exp = yper_exp/true_spectrum[1:(T/2-1)]

      # Periodogram distribution of own simulation for FARIMA
      ar_lambda[sim,j] <- fitdistr(I_ar, "exponential")[["estimate"]]
      ar_exp_1[sim,j] = ks.test(I_ar, "pexp", 1)[["p.value"]]

      exp_lambda[sim,j] <- fitdistr(I_exp, "exponential")[["estimate"]]
      exp_exp_1[sim,j] = ks.test(I_exp, "pexp", 1)[["p.value"]]
    }

  if (sim%%50==0) {
    time_end = Sys.time()
    print(abs(time_start-time_end))
    time_start = Sys.time()
  }
}

################################################################################
##----------------------------------------------------------------------------##
## 5. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. TIME
################################################################################

saveRDS(time_ar, file = paste0(path,"time_ar.RData"))
saveRDS(time_exp, file = paste0(path,"time_exp.RData"))

################################################################################
# 5.2. AVERAGE
################################################################################

saveRDS(ar_average, file = paste0(path,"ar_average.RData"))
saveRDS(exp_average, file = paste0(path,"exp_average.RData"))

################################################################################
# 5.3. TIME DOMAIN PARAMETER EXP
################################################################################

saveRDS(true_parameteres, file = paste0(path,"true_parameteres.RData"))
saveRDS(fit_c_k_ar, file = paste0(path,"fit_c_k_ar.RData"))
saveRDS(fit_c_k_exp, file = paste0(path,"fit_c_k_exp.RData"))

################################################################################
# 5.4. TIME DOMAIN PARAMETER EXP
################################################################################

saveRDS(true_ar, file = paste0(path,"true_ar.RData"))
saveRDS(fit_ar_ar, file = paste0(path,"fit_ar_ar.RData"))
saveRDS(fit_ar_exp, file = paste0(path,"fit_ar_exp.RData"))

################################################################################
# 5.5. FREQUENCY DOMAIN PARAMETER
################################################################################

saveRDS(ar_lambda, file = paste0(path,"ar_lambda.RData"))
saveRDS(exp_lambda, file = paste0(path,"exp_lambda.RData"))

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

saveRDS(ar_exp_1, file = paste0(path,"ar_exp_1.RData"))
saveRDS(exp_exp_1, file = paste0(path,"exp_exp_1.RData"))
