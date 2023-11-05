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

#DO NOT TOUCH
active_path = dirname(rstudioapi::getActiveDocumentContext()$path)
path = paste0(active_path,"/")
symbol = "d"
d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

# GENERAL SIZES
tite_size = 2.5
subtitle_size = 2
lab_size = 1.5
axis_text = 1.5
true_lab = 1.5

# LEGEND DETAILS
graph_1 = "fepxmst"
graph_2 = "fracdiff"
legend = c(graph_1, graph_2)
legend_size = 1.5
legend_line = 1
legend_line_style = 1.5

# DISPLAY SETTING
bottom_line = 5
left_line = 5.5
top_line = 1
right_line = 8
inset_num = -0.3

################################################################################
# 5.1. TIME (Figures 1,2 and 3)
################################################################################

# LOAD DATA LISTS
own_times_farima_list = readRDS(file = paste0(path,"own_times_farima_list.RData"))
own_times_fexp_list = readRDS(file = paste0(path,"own_times_fexp_list.RData"))
r_times_list = readRDS(file = paste0(path,"r_times_list.RData"))

# COMPUTE AVERAGES: DATA MATRICES
own_times_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_times_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_times_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_times_farima_matrix[wich_d,] = colMeans(own_times_farima_list[[wich_d]])
  own_times_fexp_matrix[wich_d,] = colMeans(own_times_fexp_list[[wich_d]])
  r_times_matrix[wich_d,] = colMeans(r_times_list[[wich_d]])
}

#-------------------------------------------------------------------------------
# 5.1.2. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(own_times_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_times_fexp_matrix) = paste0(POWER,"")
own_times_fexp_matrix = as.data.frame(own_times_fexp_matrix)

graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(own_times_fexp_matrix))
matplot(t(own_times_fexp_matrix), ylim=c(0,max(r_times_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = "Average running time", cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

#-------------------------------------------------------------------------------
# 5.1.3. fracdiff farima
#-------------------------------------------------------------------------------

rownames(r_times_matrix) = paste0(d_coef_vec,"")
colnames(r_times_matrix) = paste0(POWER,"")
r_times_matrix = as.data.frame(r_times_matrix)

graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(r_times_matrix))
matplot(t(r_times_matrix), ylim=c(0,max(r_times_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = "Average running time", cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.2. AVERAGE (Figures 3 and 4)
################################################################################

# DATA LISTS
own_average_farima_list = readRDS(file = paste0(path,"own_average_farima_list.RData"))
own_average_fexp_list = readRDS(file = paste0(path,"own_average_fexp_list.RData"))
r_average_list = readRDS(file = paste0(path,"r_average_list.RData"))

# DATA MATRICES
own_average_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_average_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_average_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))


for (wich_d in 1:length(d_coef_vec)) {
  own_average_farima_matrix[wich_d,] = colMeans((own_average_farima_list[[wich_d]])^2)
  own_average_fexp_matrix[wich_d,] = colMeans((own_average_fexp_list[[wich_d]])^2)
  r_average_matrix[wich_d,] = colMeans((r_average_list[[wich_d]])^2)
}

#-------------------------------------------------------------------------------
# 5.2.2. fexpmst fexpmst
#-------------------------------------------------------------------------------

rownames(own_average_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_average_fexp_matrix) = paste0(POWER,"")
own_average_fexp_matrix = as.data.frame(own_average_fexp_matrix)

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(own_average_fexp_matrix))
matplot(t(own_average_fexp_matrix), ylim=c(0,max(r_average_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$\\bar{\\hat{\\mu}}$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

#-------------------------------------------------------------------------------
# 5.2.3. fracdiff farima
#-------------------------------------------------------------------------------

rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")
r_average_matrix = as.data.frame(r_average_matrix)

graph_name = "Figure 4.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(r_average_matrix))
matplot(t(r_average_matrix), ylim=c(0,max(r_average_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$\\bar{\\hat{\\mu}}$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()


################################################################################
# 5.3. TIME DOMAIN PARAMETER d (Figures 5 and 6)
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

#-------------------------------------------------------------------------------
# 5.3.2. fexpmst fexpmst
#-------------------------------------------------------------------------------

rownames(own_long_param_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_long_param_fexp_matrix) = paste0(POWER,"")
own_long_param_fexp_matrix = as.data.frame(own_long_param_fexp_matrix)

graph_name = "Figure 5.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(own_long_param_fexp_matrix))
matplot(t(own_long_param_fexp_matrix), ylim=c(0,max(own_long_param_fexp_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$MSE(\\bar{\\hat{d}},d)$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

#-------------------------------------------------------------------------------
# 5.3.3. fexpmst fexpm
#-------------------------------------------------------------------------------

rownames(r_long_param_matrix) = paste0(d_coef_vec,"")
colnames(r_long_param_matrix) = paste0(POWER,"")
r_long_param_matrix = as.data.frame(r_long_param_matrix)

graph_name = "Figure 6.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(r_long_param_matrix))
matplot(t(r_long_param_matrix), ylim=c(0,max(r_long_param_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$MSE(\\bar{\\hat{d}},d)$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()


################################################################################
# 5.4. FREQUENCY DOMAIN PARAMETER (Tables: 16,17,18,19,20)
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

#-------------------------------------------------------------------------------
# 5.4.2. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(own_lambda_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_lambda_fexp_matrix) = paste0(POWER,"")
own_lambda_fexp_matrix = as.data.frame(own_lambda_fexp_matrix)

graph_name = "Figure 7.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(own_lambda_fexp_matrix))
matplot(t(own_lambda_fexp_matrix), ylim=c(0,max(own_lambda_fexp_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$MSE(\\hat{\\lambda},\\lambda)$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

#-------------------------------------------------------------------------------
# 5.4.3. fracdiff farima
#-------------------------------------------------------------------------------

rownames(r_lambda_matrix) = paste0(d_coef_vec,"")
colnames(r_lambda_matrix) = paste0(POWER,"")
r_lambda_matrix = as.data.frame(r_lambda_matrix)

graph_name = "Figure 8.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(r_lambda_matrix))
matplot(t(r_lambda_matrix), ylim=c(0,max(r_lambda_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX("$MSE(\\hat{\\lambda},\\lambda)$"), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT (Tables: 21,22,23,24,25)
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

#-------------------------------------------------------------------------------
# 5.5.2. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(own_exp_1_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_exp_1_fexp_matrix) = paste0(POWER,"")
own_exp_1_fexp_matrix = as.data.frame(own_exp_1_fexp_matrix)

graph_name = "Figure 9.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(own_exp_1_fexp_matrix))
matplot(t(own_exp_1_fexp_matrix), ylim=c(0,max(own_exp_1_fexp_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX(paste0("No rejected $H_0$, (%)")), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()

#-------------------------------------------------------------------------------
# 5.5.3. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(r_exp_1_matrix) = paste0(d_coef_vec,"")
colnames(r_exp_1_matrix) = paste0(POWER,"")
r_exp_1_matrix = as.data.frame(r_exp_1_matrix)

graph_name = "Figure 10.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)
row_colors <- rainbow(nrow(r_exp_1_matrix))
matplot(t(r_exp_1_matrix), ylim=c(0,max(r_exp_1_matrix)),type = "l", col = row_colors, lty = 1,lwd = 1.5, xlab = "T", ylab = TeX(paste0("No rejected $H_0$, (%)")), cex.lab = lab_size, cex.axis = axis_text, labels = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)
axis(1, at=1:8, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
       col = row_colors, lty = 1, title = "d")#, xjust = -, yjust = 1.2)

# SAVE THE GRAPH
dev.off()


