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
# 4.1. COEFFICIENTS
# 4.2. TIME
# 4.3. TIME SERIES AVERAGE
# 4.4. PHI COEFFICIENTS
# 4.5. THETA COEFFICIENTS
# 4.6. FREQUENCY DOMAIN PARAMETER
# 4.7. FREQUENCY DOMAIN GOODNESS OF FIT

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
## 3. SIMULATION PARAMETERS                                                   ##
##----------------------------------------------------------------------------##
################################################################################

PROCESS = "ARMA"
SUBPROCESS = "sampled"
POWER = 7:14
N_SIMULATIONS = 1000
ar_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)
ma_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/3. ",PROCESS,"/3.2. ",SUBPROCESS,"/")

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION (expected runtime in mac with m1 chip: 54.35928 mins)        ##
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

# LAMBDA
fit_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# P-VAL
p_val_own_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_r_exp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# Ensure no common roots between AR and MA
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

    #COEFFICIENTS
    ar_coef = ar_coef_vec[sim]
    ma_coef = ma_coef_vec[sim]

    # SIMULATION LENGTH
    T = 2^(POWER[j])

    # SPECTRAL DENSITY OF THE PROCESS
    true_spectrum = farima.spectrum(ar = ar_coef, ma = ma_coef, n.freq = T)

    # ARMA(1,1) OWN simulations
    start_own = Sys.time()
    y_own = sim.farima(ar = ar_coef, ma = ma_coef, T = T)
    end_own = Sys.time()
    time_own[sim,j] = end_own-start_own

    # AR EXISTING PACKAGE
    start_r = Sys.time()
    y_r = arima.sim(model = list(ar = ar_coef, ma = ma_coef), n = T)
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

    a = fit.farima(y_own, p = 1,q = 1)
    fit_own_phi[sim,j] = a[["ar"]]
    fit_own_theta[sim,j] = a[["ma"]]

    a = fit.farima(y_r, p = 1, q = 1)
    fit_r_phi[sim,j] = a[["ar"]]
    fit_r_theta[sim,j] = a[["ma"]]

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
  if ((sim %% 50) == 0) {
    print(sim)
  }
}
end = Sys.time()
total_time = end-start
print(total_time)

################################################################################
##----------------------------------------------------------------------------##
## 4. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 4.1. COEFFICIENTS
################################################################################

saveRDS(ar_coef_vec, file = paste0(path,"ar_coef_vec.RData"))
saveRDS(ma_coef_vec, file = paste0(path,"ma_coef_vec.RData"))

################################################################################
# 4.2. TIME
################################################################################

colnames(time_own) = POWER
colnames(time_r) = POWER

saveRDS(time_own, file = paste0(path,"time_own.RData"))
saveRDS(time_r, file = paste0(path,"time_R.RData"))

################################################################################
# 4.3. TIME SERIES AVERAGE
################################################################################

colnames(average_own) = POWER
colnames(average_r) = POWER

saveRDS(average_own, file = paste0(path,"average_own.RData"))
saveRDS(average_r, file = paste0(path,"average_r.RData"))

################################################################################
# 4.4. PHI COEFFICIENTS
################################################################################

colnames(fit_own_phi) = POWER
colnames(fit_r_phi) = POWER

saveRDS(fit_own_phi, file = paste0(path,"phi_1_own.RData"))
saveRDS(fit_r_phi, file = paste0(path,"phi_1_r.RData"))

################################################################################
# 4.5. THETA COEFFICIENTS
################################################################################

colnames(fit_own_theta) = POWER
colnames(fit_r_theta) = POWER

saveRDS(fit_own_theta, file = paste0(path,"theta_1_own.RData"))
saveRDS(fit_r_theta, file = paste0(path,"theta_1_r.RData"))

################################################################################
# 4.6. FREQUENCY DOMAIN PARAMETER
################################################################################

colnames(fit_own_exp) = POWER
colnames(fit_r_exp) = POWER

saveRDS(fit_own_exp, file = paste0(path,"lambda_own.RData"))
saveRDS(fit_r_exp, file = paste0(path,"lambda_r.RData"))

################################################################################
# 4.7. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

colnames(p_val_own_exp) = POWER
colnames(p_val_r_exp) = POWER

saveRDS(p_val_own_exp, file = paste0(path,"p.val_own.RData"))
saveRDS(p_val_r_exp, file = paste0(path,"p.val_r.RData"))
