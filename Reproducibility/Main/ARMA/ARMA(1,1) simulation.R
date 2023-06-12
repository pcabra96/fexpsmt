################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)

################################################################################
# SEED
################################################################################

set.seed(0)

################################################################################
# PARAMETERS
################################################################################

POWER = 7:14
N_SIMULATIONS = 1000
ar_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)
ma_coef_vec = runif(N_SIMULATIONS, min = -0.9, max = 0.9)

names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# SIMULATION
################################################################################

fit_own_phi_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_own_theta_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_phi_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_theta_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

fit_own_exp_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
fit_r_exp_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

p_val_own_exp_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
p_val_r_exp_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

time_own_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))
time_r_list = matrix(0,nrow = N_SIMULATIONS, ncol = length(POWER))

for (sim in 1:N_SIMULATIONS) {
  for (pow in 1:length(POWER)) {
    T = 2^(POWER[j])
    ar_coef = ar_coef_vec[sim]
    ma_coef = ma_coef_vec[sim]
    true_spectrum = farima.spectrum(ar = ,ma= , n.freq = T)

    # ARMA(1,1) OWN simulations
    start_own = Sys.time()
    y_own = sim.farima(ar = ar_coef, ma = , T = T)
    end_own = Sys.time()
    time_own[sim,j] = end_own-start_own
  }

}
