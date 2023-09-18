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
this_path = paste0(active_path,"/")
################################################################################
# PARAMETERS
################################################################################

set.seed(0)

# DISPLAY SETTING
bottom_line = 5
left_line = 7
top_line = 3
right_line = 2
inset_num = -0.45

# GENERAL SIZES
tite_size = 2.5
subtitle_size = 2
lab_size = 1.5
axis_text = 1.5
true_lab = 1.5

# PARAMAETERS
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

for (ttt in 1:length(max_vec)) {

  max = max_vec[length(max_vec)+1-ttt]
  min = min_vec[ttt]
  # Zero Matrix

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


  diff_matrix = matrix(0,nrow = N_SIMULATION, ncol = order)

  for (row in 1:N_SIMULATION) {
    coef_ar = mat[row,]
    f_t_1 = farima.spectrum(ar = coef_ar, n.freq = n)
    coef_c_k = fourier.series(log(f_t_1), k = order)[["coef"]]

    diff = rep(0,order)
    for (i in 1:order) {
      f_t_2 = fexp.spectrum(coef_c_k[1:i], n.freq = n)
      aux_diff = abs(f_t_1[1:length(w)]-f_t_2[1:length(w)])
      diff[i] = trapz(w, aux_diff)
    }

    diff_matrix[row,] = diff
  }

  graph_name = paste0("AR(10,roots from (",min,",",max,")) approx by EXP(p).pdf")
  pdf(file = paste0(this_path,graph_name), width = 8, height = 5) # The height of the plot in inches
  max_num = max(diff_matrix)
  par(mfrow=c(1,1), mar=c(bottom_line, left_line, top_line,right_line))
  main = TeX(paste0("AR(",order,") with sampled ",min," < $x_i$ <",max))
  matplot(x = 1:10, y = t(diff_matrix[,1:10]), main =main, ylim=c(0,max_num),
          type = "l", ylab = TeX(paste0("$\\int_{0}^{\\pi}|f_{AR(",order,")}(\\omega)-f_{EXP(p)}(\\omega)|d\\omega$")),
          xlab = "p", cex.lab = lab_size, cex.axis = axis_text)
  dev.off()
}



