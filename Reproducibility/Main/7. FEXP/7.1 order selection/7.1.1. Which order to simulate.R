################################################################################
##----------------------------------------------------------------------------##
## INDEX                                                                      ##
##----------------------------------------------------------------------------##
################################################################################

# 1. PACKAGES
# 2. SEED

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(latex2exp)
library(knitr)
library(pracma)
library(fracdiff)
library(devtools)
#devtools::install_github("pcabra96/fexpsmt", force = TRUE)
library(fexpsmt)
library(fitdistrplus)

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

active_path = dirname(rstudioapi::getActiveDocumentContext()$path)
path = paste0(active_path,"/")
N_SIM = 1000

the_K = 10
T = 2^14

d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)

diff_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIM, ncol = 10), simplify = FALSE)
diff_exp_ar_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIM, ncol = 10), simplify = FALSE)
diff_exp_d_list = replicate(length(d_coef_vec), matrix(0,nrow = N_SIM, ncol = 10), simplify = FALSE)

for (k in 1:length(d_coef_vec)) {
  diff = matrix(0,nrow = N_SIM, ncol = 10)
  diff_exp_ar = matrix(0,nrow = N_SIM, ncol = 10)
  diff_exp_d = matrix(0,nrow = N_SIM, ncol = 10)
  d_coef = d_coef_vec[k]
  for (j in 1:10) {
      the_vec = 1:10
      K = the_vec[j]
      for (i in 1:N_SIM) {
            ar_coef = runif(1, min = -0.7, max = 0.7)
            f_t = farima.spectrum(ar = ar_coef, d = 0, n.freq = T)
            coef_ck = fourier.series(f_t = log(f_t), k = K)$coef
            true_spectrum = fexp.spectrum(ck = coef_ck, d = d_coef,n.freq = T)
            true_spectrum[1] = 2*true_spectrum[2]
            true_spectrum[length(true_spectrum)] = 2*true_spectrum[length(true_spectrum)-1]

            mhalfm <- (T-1) %/% 2L
            w <- 2*pi/T * (1:mhalfm)
            diff[i,j] = trapz(w, abs(f_t-true_spectrum)[1:length(w)])
            if(diff[i]>20){
              print(paste0("The diff is: ",diff[i],". d = ",d_coef,", ", "ar = ",ar_coef))
            }
            ##################################
            # FITTING
            ##################################
            y_ar = sim.farima(ar = ar_coef, d = d_coef, T = T)
            y_exp = sim.fexp(ck = coef_ck, d = d_coef, T = T)

            ##################################
            # FITTING
            ##################################
            a_ar = fracdiff(x = y_ar,nar = 1)
            a_exp = fracdiff(x = y_exp,nar = 1)

            ##################################
            # COEFFICIENTS SHORT MEMORY
            ##################################
            diff_exp_ar[i,j] = abs(a_exp$ar-ar_coef)

            ##################################
            # COEFFICIENTS LONG MEMORY
            ##################################
            diff_exp_d[i,j] = abs(a_exp$d-d_coef)
        }
  }
  diff_list[[k]] = diff
  diff_exp_ar_list[[k]] = diff_exp_ar
  diff_exp_d_list[[k]] = diff_exp_d
}

################################################################################
##----------------------------------------------------------------------------##
## 5. RESULTS                                                                 ##
##----------------------------------------------------------------------------##
################################################################################

saveRDS(diff_list, file = paste0(path,"diff_list.RData"))
saveRDS(diff_exp_ar_list, file = paste0(path,"diff_exp_ar_list.RData"))
saveRDS(diff_exp_d_list, file = paste0(path,"diff_exp_d_list.RData"))
