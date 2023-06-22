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

# DATA LISTS
own_times_farima_list = readRDS(file = paste0(path,"own_times_farima_list.RData"))
own_times_fexp_list = readRDS(file = paste0(path,"own_times_fexp_list.RData"))
r_times_list = readRDS(file = paste0(path,"r_times_list.RData"))

# DATA MATRICES
own_times_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_times_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_times_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_times_farima_matrix[wich_d,] = colMeans(own_times_farima_list[[wich_d]])
  own_times_fexp_matrix[wich_d,] = colMeans(own_times_fexp_list[[wich_d]])
  r_times_matrix[wich_d,] = colMeans(r_times_list[[wich_d]])
}

r_times_matrix>own_times_fexp_matrix

# AXIS NAMES FOR LATEX
rownames(own_times_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_times_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_times_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_times_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(own_times_farima_matrix<=r_times_matrix)

# LATEX OUTPUT
own_times_farima_matrix = xtable(own_times_farima_matrix, digits = 5)
own_times_fexp_matrix = xtable(own_times_fexp_matrix, digits = 5)
r_times_matrix = xtable(r_times_matrix, digits = 5)

print(own_times_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_times_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_times_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

# DATA
own_long_param_farima_list = readRDS(file = paste0(path,"own_long_param_farima_list.RData"))
own_long_param_fexp_list = readRDS(file = paste0(path,"own_long_param_fexp_list.RData"))
r_long_param_list = readRDS(file = paste0(path,"r_long_param_list.RData"))

# DATA MATRICES
own_long_param_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_long_param_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_long_param_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_long_param_farima_matrix[wich_d,] = colMeans( (own_long_param_farima_list[[wich_d]]-d_coef_vec[wich_d])^2 )
  own_long_param_fexp_matrix[wich_d,] = colMeans((own_long_param_fexp_list[[wich_d]]-d_coef_vec[wich_d])^2)
  r_long_param_matrix[wich_d,] = colMeans((r_long_param_list[[wich_d]]-d_coef_vec[wich_d])^2)
}

# AXIS NAMES FOR LATEX
rownames(own_long_param_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_long_param_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_long_param_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_long_param_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_long_param_farima_matrix<=r_long_param_matrix)
comparison_matrix_1 = as.matrix(own_long_param_fexp_matrix<=r_long_param_matrix)

# LATEX OUTPUT
own_long_param_farima_matrix = xtable(own_long_param_farima_matrix, digits = 5)
own_long_param_fexp_matrix = xtable(own_long_param_fexp_matrix, digits = 5)
r_long_param_matrix = xtable(r_long_param_matrix, digits = 5)

print(own_long_param_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_long_param_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_long_param_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# DATA
own_lambda_farima_list = readRDS(file = paste0(path,"own_lambda_farima_list.RData"))
own_lambda_fexp_list = readRDS(file = paste0(path,"own_lambda_fexp_list.RData"))
r_lambda_list = readRDS(file = paste0(path,"r_lambda_list.RData"))

# DATA MATRICES
own_lambda_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_lambda_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_lambda_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_lambda_farima_matrix[wich_d,] = colMeans( (own_lambda_farima_list[[wich_d]]-1)^2 )
  own_lambda_fexp_matrix[wich_d,] = colMeans((own_lambda_fexp_list[[wich_d]]-1)^2)
  r_lambda_matrix[wich_d,] = colMeans((r_lambda_list[[wich_d]]-1)^2)
}

# AXIS NAMES FOR LATEX
rownames(own_lambda_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_lambda_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_lambda_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_lambda_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_lambda_farima_matrix<=r_lambda_matrix)
comparison_matrix_2 = as.matrix(own_lambda_fexp_matrix<=r_lambda_matrix)

# LATEX OUTPUT
own_lambda_farima_matrix = xtable(own_lambda_farima_matrix, digits = 5)
own_lambda_fexp_matrix = xtable(own_lambda_fexp_matrix, digits = 5)
r_lambda_matrix = xtable(r_lambda_matrix, digits = 5)

print(own_lambda_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_lambda_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_lambda_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# DATA
own_exp_1_farima_list = readRDS(file = paste0(path,"own_exp_1_farima_list.RData"))
own_exp_1_fexp_list = readRDS(file = paste0(path,"own_exp_1_fexp_list.RData"))
r_exp_1_list = readRDS(file = paste0(path,"r_exp_1_list.RData"))

# DATA MATRICES
own_exp_1_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_exp_1_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_exp_1_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (i in 1:length(d_coef_vec)) {
  own_exp_1_farima_matrix[i,] = colMeans(own_exp_1_farima_list[[i]]>=0.05)
  own_exp_1_fexp_matrix[i,] = colMeans(own_exp_1_fexp_list[[i]]>=0.05)
  r_exp_1_matrix[i,] = colMeans(r_exp_1_list[[i]]>=0.05)
}

# AXIS NAMES FOR LATEX
rownames(own_exp_1_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_exp_1_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_exp_1_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_exp_1_farima_matrix>=r_exp_1_matrix)
comparison_matrix_2 = as.matrix(own_exp_1_fexp_matrix>=r_exp_1_matrix)

plot(comparison_matrix_1)
plot(comparison_matrix_2)

# LATEX OUTPUT
own_exp_1_farima_matrix = xtable(own_exp_1_farima_matrix, digits = 5)
own_exp_1_fexp_matrix = xtable(own_exp_1_fexp_matrix, digits = 5)
r_exp_1_matrix = xtable(r_exp_1_matrix, digits = 5)

print(own_exp_1_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_exp_1_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_exp_1_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})
