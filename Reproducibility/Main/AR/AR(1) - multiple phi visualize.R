################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)
library(xtable)
library(Hmisc)
library(kableExtra)

################################################################################
# PARAMETERS
################################################################################

ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# Load data
################################################################################

path = "~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/"

fit_own_coef_list = readRDS(file = paste0(path,"phi_1_own.RData"))
fit_r_coef_list = readRDS(file = paste0(path,"phi_1_r.RData"))
fit_own_exp_list = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp_list = readRDS(file = paste0(path,"lambda_r.RData"))
p_val_own_exp_list = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp_list = readRDS(file = paste0(path,"p.val:R.RData"))
time_own_list = readRDS(file = paste0(path,"time_own.RData"))
time_r_list = readRDS(file = paste0(path,"\time_R.RData"))

################################################################################
# Running time
################################################################################

average_time_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
average_time_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  average_time_own[i,] = colMeans((time_own_list[[i]]))
  average_time_r[i,] = colMeans((time_r_list[[i]]))
}
rownames(average_time_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(average_time_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(average_time_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(average_time_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(average_time_own<=average_time_r)

################################################################################
# MSE COEFFICIENTS for phi
################################################################################

MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_coef_list[[i]]-ar_coef_vec[i])^2)
  MSE_r[i,] = colMeans((fit_r_coef_list[[i]]-ar_coef_vec[i])^2)
}

rownames(MSE_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(MSE_own<=MSE_r)

MSE_own = xtable(MSE_own, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)

print(MSE_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# OWN CODE
MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_exp_list[[i]]-1)^2)
  MSE_r[i,] = colMeans((fit_r_exp_list[[i]]-1)^2)
}

rownames(MSE_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(MSE_own<=MSE_r)

MSE_own = xtable(MSE_own, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)

print(MSE_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# OWN CODE
prop_non_rejection_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
prop_non_rejection_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  prop_non_rejection_own[i,] = colMeans(p_val_own_exp_list[[i]]>=0.05)
  prop_non_rejection_r[i,] = colMeans(p_val_r_exp_list[[i]]>=0.05)
}

rownames(prop_non_rejection_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(prop_non_rejection_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(prop_non_rejection_own<prop_non_rejection_r)

prop_non_rejection_own = xtable(prop_non_rejection_own, digits = 5)
prop_non_rejection_r = xtable(prop_non_rejection_r, digits = 5)

print(prop_non_rejection_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(prop_non_rejection_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})
