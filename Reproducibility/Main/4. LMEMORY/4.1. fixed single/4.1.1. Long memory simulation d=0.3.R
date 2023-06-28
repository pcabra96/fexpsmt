################################################################################
# INDEX
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
# 1. PACKAGES
################################################################################

library(fracdiff)
library(latex2exp)

################################################################################
# 2. SEED
################################################################################

set.seed(0)

################################################################################
# 3. SIMULATION PARAMETERS
################################################################################

PROCESS = "LMEMORY"
SUBPROCESS = "fixed"
d_coef = 0.3
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
symbol = "d"

################################################################################
# Save data
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/4. ",PROCESS,"/",SUBPROCESS,"/")

# TIME
own_times_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
own_times_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
r_times = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

# AVERAGE
own_average_farima = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
own_average_fexp = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
r_average = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

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

# F_ARMA(0,0)
f_0 = farima.spectrum(ar = 0, ma = 0, n.freq = 1)

################################################################################
# 4. SIMULATION
################################################################################

start = Sys.time()

for (sim in 1:N_SIMULATIONS) {
  for (j in 1:length(POWER)) {

    T = 2^POWER[j]

    true_spectrum = farima.spectrum(d = d_coef, n.freq = T)
    frequency = c(0,pi)
    freq_2 = seq(frequency[1],frequency[2],length.out = T+1)
    lamda_aprox = (freq_2[2]-freq_2[1])/T
    true_spectrum[1] = f_0*abs(lamda_aprox)^(-2*d_coef)

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
    y_R = fracdiff.sim(n = as.numeric(T), d = d_coef)[["series"]]
    r_end = Sys.time()
    r_times[sim,j] = r_end-r_start

    ################################################################################
    # AVERAGE
    ################################################################################

    own_average_farima[sim,j] = mean(y_own_farima)
    own_average_fexp[sim,j] = mean(y_own_fexp)
    r_average[sim,j] = mean(y_R)

    ##################################
    # FITTING
    ##################################

    own_long_param_farima[sim,j] = fracdiff(x = y_own_farima, nar = 0, nma = 0)[["d"]]
    own_long_param_fexp[sim,j] = fracdiff(x = y_own_fexp, nar = 0, nma = 0)[["d"]]
    r_long_param[sim,j] = fracdiff(x = y_R, nar = 0, nma = 0)[["d"]]

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
end = Sys.time()
print(end-start)

################################################################################
# Running time
################################################################################

colnames(own_times_farima) = POWER
colnames(own_times_fexp) = POWER
colnames(r_times) = POWER

saveRDS(own_times_farima, file = paste0(path,"own_times_farima.RData"))
saveRDS(own_times_fexp, file = paste0(path,"own_times_fexp.RData"))
saveRDS(r_times, file = paste0(path,"r_times.RData"))

################################################################################
# Time Series average
################################################################################

colnames(own_average_farima) = POWER
colnames(own_average_fexp) = POWER
colnames(r_average) = POWER

saveRDS(own_average_farima, file = paste0(path,"own_average_farima.RData"))
saveRDS(own_average_fexp, file = paste0(path,"own_average_fexp.RData"))
saveRDS(r_average, file = paste0(path,"r_average.RData"))

################################################################################
# MSE COEFFICIENTS for d
################################################################################

colnames(own_long_param_farima) = POWER
colnames(own_long_param_fexp) = POWER
colnames(r_long_param) = POWER

saveRDS(own_long_param_farima, file = paste0(path,"own_long_param_farima.RData"))
saveRDS(own_long_param_fexp, file = paste0(path,"own_long_param_fexp.RData"))
saveRDS(r_long_param, file = paste0(path,"r_long_param.RData"))

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
colnames(own_lambda_farima) = POWER
colnames(own_lambda_fexp) = POWER
colnames(r_lambda) = POWER

saveRDS(own_lambda_farima, file = paste0(path,"own_lambda_farima.RData"))
saveRDS(own_lambda_fexp, file = paste0(path,"own_lambda_fexp.RData"))
saveRDS(r_lambda, file = paste0(path,"r_lambda.RData"))

################################################################################
# Correct p-values for EXP(1)
################################################################################

# OWN CODE
colnames(own_exp_1_farima) = POWER
colnames(own_exp_1_fexp) = POWER
colnames(r_exp_1) = POWER

saveRDS(own_exp_1_farima, file = paste0(path,"own_exp_1_farima.RData"))
saveRDS(own_exp_1_fexp, file = paste0(path,"own_exp_1_fexp.RData"))
saveRDS(r_exp_1, file = paste0(path,"r_exp_1.RData"))
