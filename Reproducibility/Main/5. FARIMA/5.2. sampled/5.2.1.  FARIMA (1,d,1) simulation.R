################################################################################
# INDEX
################################################################################

# 1. PACKAGES
# 2. SEED
# 3. SIMULATION PARAMETERS
# 4. SIMULATION
# 5. RESULTS
# 5.1. COEFFICIENTS
# 5.2. TIME
# 5.3. TIME SERIES AVERAGE
# 5.4. PHI COEFFICIENTS
# 5.5. THETA COEFFICIENTS
# 5.6. d COEFFICIENTS
# 5.7. FREQUENCY DOMAIN PARAMETER
# 5.8. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(fracdiff)
library(forecast)
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

PROCESS = "FARIMA"
SUBPROCESS = "sampled"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/5. ",PROCESS,"/5.2. ",SUBPROCESS,"/")
POWER = 7:14
N_SIMULATIONS = 1000
ar_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)
ma_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)
d_coef_vec = runif(N_SIMULATIONS, min = 0, max = 0.5)

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION (expected runtime in mac with m1 chip:  9.8 mins)            ##
##----------------------------------------------------------------------------##
################################################################################

# TIME
time_own = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
time_r = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# AVERAGE
average_own = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
average_r = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# AR COEF
fit_own_phi = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_own_theta = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# MA COEF
fit_r_phi = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_theta = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# d coefficient
fit_r_d = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_own_d = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# LAMBDA
fit_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# P-VAL
p_val_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

for (sim in 1:N_SIMULATIONS) {
  if(abs(ma_coef_vec[sim]+ar_coef_vec[sim])<0.1){
    if (abs(ar_coef_vec[sim])<0.1 ){
      ma_coef_vec[sim] = runif(1, min = 0.2, max = 0.9)
    }
    else{
      ma_coef_vec[sim] = ar_coef_vec[sim]
    }
  }
}

start = Sys.time()

for (sim in 1:N_SIMULATIONS) {
  for (j in 1:length(POWER)) {
    T = 2^(POWER[j])

    ar_coef = ar_coef_vec[sim]
    ma_coef = ma_coef_vec[sim]
    d_coef = d_coef_vec[sim]
    true_spectrum = farima.spectrum(ar = ar_coef,ma= ma_coef, d = d_coef, n.freq = T)
    true_spectrum[1] = true_spectrum[length(true_spectrum)]

    # FARIMA(1,d,1) OWN simulations
    start_own = Sys.time()
    y_own = sim.farima(ar = ar_coef, ma = ma_coef, d = d_coef, T = T)
    end_own = Sys.time()
    time_own[sim,j] = end_own-start_own

    # FARIMA(1,d,1) fracdiff simulations
    start_r = Sys.time()
    y_r = fracdiff.sim(n = T, ar = ar_coef, ma = -ma_coef, d = d_coef)[["series"]]
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

    a = fracdiff(x = y_own,nar = 1, nma = 1)
    fit_own_phi[sim,j] = a[["ar"]]
    fit_own_theta[sim,j] = -a[["ma"]]
    fit_own_d[sim,j] = a[["d"]]

    a = fracdiff(x = y_r,nar = 1, nma = 1)
    fit_r_phi[sim,j] = a[["ar"]]
    fit_r_theta[sim,j] = -a[["ma"]]
    fit_r_d[sim,j] = a[["d"]]

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

  if ((sim %% 10) == 0) {
    print(sim)
  }
}

end = Sys.time()
total_time = end-start
print(total_time)

################################################################################
##----------------------------------------------------------------------------##
## 5. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. COEFFICIENTS
################################################################################

saveRDS(ar_coef_vec, file = paste0(path,"ar_coef_vec.RData"))
saveRDS(ma_coef_vec, file = paste0(path,"ma_coef_vec.RData"))
saveRDS(d_coef_vec, file = paste0(path,"d_coef_vec.RData"))

################################################################################
# 5.2. TIME
################################################################################

colnames(time_own) = POWER
colnames(time_r) = POWER

saveRDS(time_own, file = paste0(path,"time_own.RData"))
saveRDS(time_r, file = paste0(path,"time_R.RData"))

################################################################################
# 5.3. TIME SERIES AVERAGE
################################################################################

colnames(average_own) = POWER
colnames(average_r) = POWER

saveRDS(average_own, file = paste0(path,"average_own.RData"))
saveRDS(average_r, file = paste0(path,"average_r.RData"))

################################################################################
# 5.4. PHI COEFFICIENTS
################################################################################

colnames(fit_own_phi) = POWER
colnames(fit_r_phi) = POWER

saveRDS(fit_own_phi, file = paste0(path,"phi_1_own.RData"))
saveRDS(fit_r_phi, file = paste0(path,"phi_1_r.RData"))

################################################################################
# 5.5. THETA COEFFICIENTS
################################################################################

colnames(fit_own_theta) = POWER
colnames(fit_r_theta) = POWER

saveRDS(fit_own_theta, file = paste0(path,"theta_1_own.RData"))
saveRDS(fit_r_theta, file = paste0(path,"theta_1_r.RData"))

################################################################################
# 5.6. d COEFFICIENTS
################################################################################

colnames(fit_own_d) = POWER
colnames(fit_r_d) = POWER

saveRDS(fit_own_d, file = paste0(path,"d_own.RData"))
saveRDS(fit_r_d, file = paste0(path,"d_r.RData"))

################################################################################
# 5.7. FREQUENCY DOMAIN PARAMETER
################################################################################

# OWN CODE
colnames(fit_own_exp) = POWER
colnames(fit_r_exp) = POWER

saveRDS(fit_own_exp, file = paste0(path,"lambda_own.RData"))
saveRDS(fit_r_exp, file = paste0(path,"lambda_r.RData"))

################################################################################
# 5.8. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

# OWN CODE
colnames(p_val_own_exp) = POWER
colnames(p_val_r_exp) = POWER

saveRDS(p_val_own_exp, file = paste0(path,"p.val_own.RData"))
saveRDS(p_val_r_exp, file = paste0(path,"p.val_r.RData"))
