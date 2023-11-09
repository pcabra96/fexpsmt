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

library(devtools)
devtools::install_github("pcabra96/fexpsmt", force = TRUE)
library(fexpsmt)

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

# GENERAL SIZES
tite_size = 2.5
subtitle_size = 2
lab_size = 2
axis_text = 2
true_lab = 2

# LEGEND DETAILS
legend_size = 3
legend_line = 1
legend_line_style = 1.75

# DISPLAY SETTING
bottom_line = 5
left_line = 7.5
top_line = 1
right_line = 8
inset_num = -0.3

# PARAMETERS
N_SIM = 100
T = 2^12
mhalfm <- (T-1) %/% 2L
w <- 2*pi/T * (1:mhalfm)
k = 1:10

################################################################################
# phi1 between -0.9, 0.9
################################################################################

vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

diff = matrix(0,nrow = length(vec), ncol = 10)
for (j in 1:10) {
  for (i in 1:length(vec)) {
      ar_spectrum = farima.spectrum(ar = vec[i],n.freq = T)
      exp_coef = fourier.series(log(ar_spectrum), k = j)$coef
      exp_spectrum = fexp.spectrum(ck = exp_coef, n.freq = T)
      aux_diff = abs(ar_spectrum[1:length(w)]-exp_spectrum[1:length(w)])
      diff[i,j] = trapz(w, aux_diff)
  }
}

row.names(diff) = vec
colnames(diff) = 1:10


graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(diff))
matplot(t(diff), ylim=c(0,max(diff)), type = "l", col = row_colors, xlab = "p'", ylab = TeX("$\\int_{0}^{\\pi}|f_{AR(\\phi_1)}(\\omega)-f_{EXP(p)}(\\omega)|\\partial \\omega$"),
        lty = 1,lwd = 1.5, cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
axis(1, at=1:10, labels = 1:10)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = vec, col = row_colors, lty = 1, title = TeX("$\\phi_1$"))

# SAVE THE GRAPH
dev.off()



vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

diff_coef = matrix(0,nrow = length(vec), ncol = 10)
for (j in 1:10) {
  for (i in 1:length(vec)) {
    ar_spectrum = farima.spectrum(ar = vec[i],n.freq = T*4)
    exp_coef = fourier.series(log(ar_spectrum), k = j)$coef
    y_2 = sim.fexp(exp_coef, T = T*4)
    ar_1_hat = fit.farima(y_2, p = 1)$ar
    diff_coef[i,j] = abs(vec[i]-ar_1_hat)
  }
}



row.names(diff_coef) = vec
colnames(diff_coef) = 1:10


graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(diff_coef))
matplot(t(diff_coef), ylim=c(0,max(diff_coef)), type = "l", col = row_colors, xlab = "p'", ylab = TeX("$|\\hat{\\phi}_{1,Whittle}-\\phi_1|$"),
        lty = 1,lwd = 1.5, cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
axis(1, at=1:10, labels = 1:10)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = vec, col = row_colors, lty = 1, title = TeX("$\\phi_1$"))

# SAVE THE GRAPH
dev.off()


