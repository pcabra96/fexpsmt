################################################################################
# PACKAGES
################################################################################

library(latex2exp)
library(polynom)
library(pracma)
library(RColorBrewer)
library(gt)
library(gtExtras)
library(magrittr)
library(htmltools)
library(ggplot2)
library(tidyr)
library(reshape2)
library(cowplot)
library(ggnewscale)

active_path = dirname(rstudioapi::getActiveDocumentContext()$path)

################################################################################
# PARAMETERS
################################################################################

set.seed(0)

# GENERAL
tite_size = 2.5
subtitle_size = 2
lab_size = 1.8
axis_text = 1.5
true_lab = 1.5
bottom_box = 5
left_box = 6
top_box = 6
right_box = 2

n = 2^12
mhalfm <- (n-1) %/% 2L
w <- 2*pi/n * (1:mhalfm)

N_SIMULATION = 250
order = 10
min_vec = seq(-0.9,-0.1,0.1)
max_vec = seq(0.1,0.9,0.1)

################################################################################
# RANDOM POLYNOMIALS COEFFICIENTS
################################################################################


  max = 0.3
  min = -0.3

  # Zero Matrix
  #if (max == 0 & min<0) {
  #  this_path = paste0(active_path, "/Negative roots/")
  #}
  #if (max>0 & min == 0) {
  #this_path = paste0(active_path, "/Positive roots/")
  #}
  #if (max > 0 & min < 0) {
  #    this_path = paste0(active_path, "/Mixed roots/")
  #}

  mat = matrix(0,nrow = N_SIMULATION, ncol = order)

  # Loop to store the coefficients of each AR(10)
  for (j in 1:N_SIMULATION) {
    uniform_sample = runif(order,min = min, max = max)
    # Loop to get the 10 different coefficients
    for (pol_order in 1:(length(uniform_sample))) {
      if (pol_order==1) {
        pr=1
      }
      pr = pr*polynomial(c(1,uniform_sample[pol_order]))
    }
    #Store the polynomial coefficients
    pol_coef = -coefficients(pr)[-1]
    mat[j,] = pol_coef
  }


  ar_coef = mat[1,]
  y = sim.farima(ar = ar_coef, T = 2^13)
  f_t = farima.spectrum(ar = ar_coef, n.freq = 2^13)
  coef_c_k = fourier.series(f_t = log(f_t), k=3)[["coef"]]
  y_2 = sim.fexp(ck = coef_c_k, T = 2^13)

  par(mfrow=c(2,1))
  plot(y[1:512], main = "AR(10)", ylab = "", type = "l")
  plot(y_2[1:512], main = "EXP(3)", ylab = "", type = "l")

  par(mfrow=c(2,1))
  acf(y, main = "AR(10)")
  acf(y_2, main ="EXP(3)")

  fit_ar_ar = fit.farima(y, p=10)$ar
  fit_ar_exp = fit.farima(y_2, p=10)$ar

  fit_exp_ar = fit.fexp(y, p=2)$c_k
  fit_exp_exp = fit.fexp(y_2, p=2)$c_k

