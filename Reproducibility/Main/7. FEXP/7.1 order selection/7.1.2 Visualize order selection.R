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
# 5.1. TIME (Figures 1,2,3)
# 5.2. AVERAGE (Figures 2,3,4)
# 5.3. TIME DOMAIN PARAMETER d (Figures 4,5,6)
# 5.4. FREQUENCY DOMAIN PARAMETER (7,8,9)
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT (10,11,12)

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(RColorBrewer)
library(latex2exp)
library(gt)
library(gtExtras)
library(magrittr)
library(htmltools)
library(ggplot2)
library(tidyr)
library(reshape2)

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

# GENERAL
tite_size = 2.5
subtitle_size = 2
lab_size = 1.5
axis_text = 1.5
true_lab = 1.5

# DISPLAY SETTING
bottom_line = 5
left_line = 7
top_line = 1
right_line = 7.5
inset_num = -0.25

d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)

diff_list = readRDS(file = paste0(path,"diff_list.RData"))
diff_exp_ar_list = readRDS(file = paste0(path,"diff_exp_ar_list.RData"))
diff_exp_d_list = readRDS(file = paste0(path,"diff_exp_d_list.RData"))

################################################################################
# Spectral density difference
################################################################################

diff_to_plot = matrix(0,nrow = length(diff_list), ncol = ncol(diff_list[[1]]))
for (i in 1:length(diff_list)) {
  print(i)
  diff_to_plot[i,] = colMeans(diff_list[[i]])
}

diff_to_plot = t(diff_to_plot)
colnames(diff_to_plot) <- d_coef_vec

graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(length(d_coef_vec))
matplot(diff_to_plot, ylim=c(0,max(diff_to_plot)), type = "l", xlab = "p'", col = row_colors,
        ylab  = TeX("$\\int_{0}^{\\pi}|f_{FEXP(p,d)}(\\omega)-f_{FARIMA(1,d,0)}(\\omega)|\\partial\\omega$"), lty = 1,lwd = 1.5, cex.lab = lab_size, cex.axis = axis_text) #, type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "p'", ylab = TeX("$|\\hat{\\phi}_{1,Whittle}-\\phi_1|$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
legend("topright", inset = c(inset_num, 0), legend = d_coef_vec, col = row_colors, lty = 1, title = "d")

# SAVE THE GRAPH
dev.off()

################################################################################
# AR Difference
################################################################################

diff_exp_ar_to_plot = matrix(0,nrow = length(diff_list), ncol = ncol(diff_list[[1]]))
for (i in 1:length(diff_list)) {
  diff_exp_ar_to_plot[i,] = colMeans(diff_exp_ar_list[[i]])
}

diff_exp_ar_to_plot = t(diff_exp_ar_to_plot)
colnames(diff_exp_ar_to_plot) <- d_coef_vec

graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(length(d_coef_vec))
matplot(diff_exp_ar_to_plot, ylim=c(0,max(diff_exp_ar_to_plot)), type = "l", xlab = "p'", col = row_colors,
        ylab  = TeX("$|\\hat{\\phi}_1-\\phi_1|$"), lty = 1,lwd = 1.5, cex.lab = lab_size, cex.axis = axis_text)
legend("topright", inset = c(inset_num, 0), legend = d_coef_vec, col = row_colors, lty = 1, title = "d",
       cex = legend_size,lwd=legend_line_style,seg.len=legend_line)

# SAVE THE GRAPH
dev.off()

################################################################################
# d DIFFERENCE
################################################################################

diff_exp_d_to_plot = matrix(0,nrow = length(diff_list), ncol = ncol(diff_list[[1]]))
for (i in 1:length(diff_list)) {
  diff_exp_d_to_plot[i,] = colMeans(diff_exp_d_list[[i]])
}

diff_exp_d_to_plot = t(diff_exp_d_to_plot)
colnames(diff_exp_d_to_plot) <- d_coef_vec

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(length(d_coef_vec))
matplot(diff_exp_d_to_plot, ylim=c(0,max(diff_exp_d_to_plot)), type = "l", xlab = "p'", col = row_colors,
        ylab  = TeX("$|\\hat{d}-d|$"), lty = 1,lwd = 1.5, cex.lab = lab_size, cex.axis = axis_text)

legend("topright", inset = c(inset_num, 0), legend = d_coef_vec, col = row_colors, lty = 1, title = "d",
       cex = legend_size,lwd=legend_line_style,seg.len=legend_line)

# SAVE THE GRAPH
dev.off()
