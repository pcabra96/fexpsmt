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
# 5.2. TIME SERIES AVERAGE
# 5.3. PHI COEFFICIENT
# 5.3. FREQUENCY DOMAIN PARAMETER
# 5.4. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(forecast)
require(MASS)
library(latex2exp)

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

PROCESS = "MA"
SUBPROCESS = "fixed single"
symbol = "\\phi"
coef = -0.7
POWER = 7:14
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/2. ",PROCESS,"/2.1. ",SUBPROCESS,"/")

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION                   3.408257 mins                                           ##
##----------------------------------------------------------------------------##
################################################################################

N_SIMULATIONS = 1000
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

time_start = Sys.time()
for (sim in 1:N_SIMULATIONS) {
  for (j in 1:length(POWER)) {
    T = 2^(POWER[j])
    true_spectrum = farima.spectrum(ma = coef, n.freq = T)

    # AR OWN simulations
    own_start = Sys.time()
    y_own = sim.farima(ma = coef, T = T)
    own_end = Sys.time()
    time_own[sim,j] = own_end - own_start

    # AR EXISTING PACKAGE
    r_start = Sys.time()
    y_r = arima.sim(model = list(ma = coef), n = T)
    r_end = Sys.time()
    time_r[sim,j] = r_end - r_start

    ################################################################################
    # AVERAGE
    ################################################################################

    average_own[sim,j] = mean(y_own)
    average_r[sim,j] = mean(y_r)
    ################################################################################
    # FITTING
    ################################################################################

    fit_own_coef[sim,j] = fit.farima(y_own, q = 1)[["ma"]]
    fit_r_coef[sim,j] = fit.farima(y_r, q = 1)[["ma"]]

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
  if (sim %% 50 ==0) {
    print(sim)
    if (sim==50) {
      aux_time = Sys.time()
      t_aux_to_print = aux_time - time_start
      print(t_aux_to_print)
    }else{
      aux_time_2 = Sys.time()
      t_aux_to_print = aux_time_2 - aux_time
      print(t_aux_to_print)
    }
  }
}
time_end = Sys.time()
time = time_end-time_start
print(time)

################################################################################
##----------------------------------------------------------------------------##
## 5. SAVE RESULTS                                                            ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. Running time
################################################################################

colnames(time_own) = POWER
colnames(time_r) = POWER

saveRDS(time_own, file = paste0(path,"time_own.RData"))
saveRDS(time_r, file = paste0(path,"time_R.RData"))

################################################################################
# 5.2. Time Series average
################################################################################

colnames(average_own) = POWER
colnames(average_r) = POWER

saveRDS(average_own, file = paste0(path,"average_own.RData"))
saveRDS(average_r, file = paste0(path,"average_r.RData"))

################################################################################
# 5.3. MSE COEFFICIENTS for phi
################################################################################

colnames(fit_own_coef) = POWER
colnames(fit_r_coef) = POWER

saveRDS(fit_own_coef, file = paste0(path,"theta_1_own.RData"))
saveRDS(fit_r_coef, file = paste0(path,"theta_1_r.RData"))

################################################################################
# 5.4. MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
colnames(fit_own_exp) = POWER
colnames(fit_r_exp) = POWER

saveRDS(fit_own_exp, file = paste0(path,"lambda_own.RData"))
saveRDS(fit_r_exp, file = paste0(path,"lambda_r.RData"))

################################################################################
# 5.5. Correct p-values for EXP(1)
################################################################################

# OWN CODE
colnames(p_val_own_exp) = POWER
colnames(p_val_r_exp) = POWER

saveRDS(p_val_own_exp, file = paste0(path,"p.val_own.RData"))
saveRDS(p_val_r_exp, file = paste0(path,"p.val_r.RData"))
