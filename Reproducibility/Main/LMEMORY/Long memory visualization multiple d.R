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

library(arfima)
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

PROCESS = "FARIMA"
symbol = "d"
d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# PATH TO LOAD THE DATA
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/LMEMORY/")

################################################################################
# Running time
################################################################################

# DATA
own_times_farima_list = readRDS(file = paste0(path,"own_times_farima_list.RData"))
own_times_fexp_list = readRDS(file = paste0(path,"own_times_fexp_list.RData"))
r_times_list = readRDS(file = paste0(path,"r_times_list.RData"))

# AXIS NAMES FOR LATEX
rownames(time_own_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(time_own_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(time_r_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(time_r_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(time_own_matrix<=time_r_matrix)

# LATEX OUTPUT
time_own_matrix = xtable(time_own_matrix, digits = 5)
time_r_matrix = xtable(time_r_matrix, digits = 5)
print(time_own_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(time_own_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(time_r_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(time_r_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

# DATA
fit_own_coef_matrix = readRDS(file = paste0(path,"theta_1_own.RData"))
fit_r_coef_matrix = readRDS(file = paste0(path,"theta_1_r.RData"))

# AXIS NAMES FOR LATEX
rownames(fit_own_coef_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(fit_own_coef_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(fit_r_coef_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(fit_r_coef_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(fit_own_coef_matrix<=fit_r_coef_matrix)

# LATEX OUTPUT
fit_own_coef_matrix = xtable(fit_own_coef_matrix, digits = 5)
fit_r_coef_matrix = xtable(fit_r_coef_matrix, digits = 5)
print(fit_own_coef_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(fit_own_coef_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(fit_r_coef_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(fit_r_coef_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# DATA
fit_own_exp_matrix = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp_matrix = readRDS(file = paste0(path,"lambda_r.RData"))

# AXIS NAMES FOR LATEX
rownames(fit_own_exp_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(fit_own_exp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(fit_r_exp_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(fit_r_exp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(fit_own_exp_matrix<=fit_r_exp_matrix)

# LATEX OUTPUT
fit_own_exp_matrix = xtable(fit_own_exp_matrix, digits = 5)
fit_r_exp_matrix = xtable(fit_r_exp_matrix, digits = 5)
print(fit_own_exp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(fit_own_exp_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(fit_r_exp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(fit_r_exp_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# DATA
p_val_own_exp_matrix = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp_matrix = readRDS(file = paste0(path,"p.val_R.RData"))

# AXIS NAMES FOR LATEX
rownames(p_val_own_exp_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(p_val_own_exp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(p_val_r_exp_matrix) = paste0("$\\mathbf{",ma_coef_vec,"}$")
colnames(p_val_r_exp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(p_val_own_exp_matrix<p_val_r_exp_matrix)

# LATEX OUTPUT
p_val_own_exp_matrix = xtable(p_val_own_exp_matrix, digits = 5)
p_val_r_exp_matrix = xtable(p_val_r_exp_matrix, digits = 5)
print(p_val_own_exp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(p_val_own_exp_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(p_val_r_exp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(p_val_r_exp_matrix)), sanitize.text.function = function(x) {x})
