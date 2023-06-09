################################################################################
# PACKAGES
################################################################################

library(forecast)
require(MASS)
library(latex2exp)
library(xtable)
library(Hmisc)
library(kableExtra)
library(plot.matrix)

################################################################################
# PARAMETERS
################################################################################

ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# PATH TO LOAD THE DATA
################################################################################

path = "~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/AR/"

################################################################################
# Running time
################################################################################

# DATA
time_own_matrix = readRDS(file = paste0(path,"time_own.RData"))
time_r_matrix = readRDS(file = paste0(path,"time_R.RData"))

# AXIS NAMES FOR LATEX
rownames(time_own_matrix) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(time_own_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(time_r_matrix) = paste0("$\\mathbf{",ar_coef_vec,"}$")
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
fit_own_coef_matrix = readRDS(file = paste0(path,"phi_1_own.RData"))
fit_r_coef_matrix = readRDS(file = paste0(path,"phi_1_r.RData"))

# AXIS NAMES FOR LATEX
rownames(fit_own_coef_matrix) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(fit_own_coef_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(fit_own_coef_matrix<=MSE_r)

# LATEX OUTPUT
fit_own_coef_matrix = xtable(fit_own_coef_matrix, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)
print(fit_own_coef_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(fit_own_coef_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# DATA
fit_own_exp_matrix = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp_matrix = readRDS(file = paste0(path,"lambda_r.RData"))

# AXIS NAMES FOR LATEX
rownames(MSE_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(MSE_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(MSE_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(MSE_own<=MSE_r)

# LATEX OUTPUT
MSE_own = xtable(MSE_own, digits = 5)
MSE_r = xtable(MSE_r, digits = 5)
print(MSE_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(MSE_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# DATA
p_val_own_exp_matrix = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp_matrix = readRDS(file = paste0(path,"p.val_R.RData"))

# AXIS NAMES FOR LATEX
rownames(prop_non_rejection_own) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_own) = paste0("$\\mathbf{2^{",POWER,"}}$")
rownames(prop_non_rejection_r) = paste0("$\\mathbf{",ar_coef_vec,"}$")
colnames(prop_non_rejection_r) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix = as.matrix(prop_non_rejection_own<prop_non_rejection_r)

# LATEX OUTPUT
prop_non_rejection_own = xtable(prop_non_rejection_own, digits = 5)
prop_non_rejection_r = xtable(prop_non_rejection_r, digits = 5)
print(prop_non_rejection_own, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_own)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(prop_non_rejection_r, include.rownames = TRUE, hline.after = c(-1,0, nrow(MSE_r)), sanitize.text.function = function(x) {x})
