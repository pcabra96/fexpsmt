################################################################################
##----------------------------------------------------------------------------##
## INDEX                                                                      ##
##----------------------------------------------------------------------------##
################################################################################

# 1. PACKAGES
# 2. SEED
# 3. SIMULATION PARAMETERS
# 4. VISUALIZATION PATH
# 5. VISUALIZATION
# 5.1. TIME (Figures: 1,2)
# 5.2. AVERAGE (Figures: 3,4)
# 5.3. MSE FOR AR COEFFICIENT (Figures: 5,6)
# 5.4. MSE FOR MA COEFFICIENT (Figures: 7,8)
# 5.5. MSE FOR d COEFFICIENT (Figures: 9,10)
# 5.6. MSE FOR LAMBDA COEFFICIENT (Figures: 11,12)
# 5.7. FREQUENCY DOMAIN GOODNESS OF FIT (Figures: 13,14)

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

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

PROCESS = "EXP"
SUBPROCESS = "sampled"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/6. ",PROCESS,"/6.3. ",SUBPROCESS,"/")
POWER = 7:14
K = 4

# DO NOT TOUCH
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14
N_SIM = 1000

################################################################################
##----------------------------------------------------------------------------##
## 4. VISUALIZATION PARAMETERS                                                ##
##----------------------------------------------------------------------------##
################################################################################

# GENERAL
tite_size = 2.5
subtitle_size = 2
lab_size = 1.8
axis_text = 1.5
true_lab = 1.5

graph_1 = TeX("$fepxmst_{AR(1)}$")
graph_2 = TeX("$fepxmst_{EXP(4)}$")
legend = c(graph_1, graph_2)
legend_size = 2
legend_line = 1
legend_line_style = 1.5

width_graphs = 1100
height_graphs = 550

# BOXPLOT c(bottom, left, top, right))mar=c(5,5,6,2)
bottom_box = 5
left_box = 6
top_box = 6
right_box = 2

# LINES c(bottom, left, top, right))mar=c(5,5,6,2)
bottom_line = 5
left_line = 6
top_line = 6
right_line = 2

################################################################################
##----------------------------------------------------------------------------##
## 5. VISUALIZATION                                                           ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. TIME
################################################################################

dev.off()
graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# LOAD THE DATA
time_ar = readRDS(file = paste0(path,"time_ar.RData"))
time_exp = readRDS(file = paste0(path,"time_exp.RData"))

#-------------------------------
# 5.1.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = ""
y_1 = time_ar
y_2 = time_exp
abline_value = 0
abline_col = ""

# LIMITS
lim_inf = min(y_1,y_2)
lim_sup = max(y_1,y_2)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,3,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

#-------------------------------
# 5.1.2. LINES
#-------------------------------

graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = ""
y_1 = colMeans(time_ar)
y_2 = colMeans(time_exp)
abline_value = 0
abline_col = ""

# LIMITS
lim_inf = min(y_1,y_2)
lim_sup = max(y_1,y_2)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topleft", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.2. Time Series average
################################################################################

# LOAD THE DATA
ar_average <- readRDS(file = paste0(path, "ar_average.RData"))
exp_average <- readRDS(file = paste0(path, "exp_average.RData"))

#-------------------------------
# 5.2.1. BOXPLOT
#-------------------------------

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{mu}}$")
xlab = "T"
main = ""
y_1 = ar_average
y_2 = exp_average
abline_value = 0
abline_col = "red"

# LIMITS
lim_inf = min(y_1,y_2)
lim_sup = max(y_1,y_2)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,3,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\mu$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\mu$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

#-------------------------------
# 5.2.2. LINES
#-------------------------------

graph_name = "Figure 4.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{mu}}$")
xlab = "T"
main = ""
true_param = 0
own_param = ar_average
r_param = exp_average
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param ))
y_2 = colMeans((r_param))

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param)), 2, sd)
sd_r <- apply(((r_param)), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.3. MSE COEFFICIENTS FOR CK COEFFICIENT
################################################################################

# LOAD THE DATA
true_parameteres <- readRDS(file = paste0(path, "true_parameteres.RData"))
fit_c_k_ar <- readRDS(file = paste0(path, "fit_c_k_ar.RData"))
fit_c_k_exp <- readRDS(file = paste0(path, "fit_c_k_exp.RData"))

# FITTING PLOTS EXP COEFFICIENTS
results_c_k_ar = matrix(0,nrow = N_SIM, ncol = length(POWER))
for (j in 1:length(POWER)) {
  results_c_k_ar[,j] = colSums(t((fit_c_k_ar[[j]]-true_parameteres[[j]])^2))
}

results_c_k_exp = matrix(0,nrow = N_SIM, ncol = length(POWER))
for (j in 1:length(POWER)) {
  results_c_k_exp[,j] = colSums(t((fit_c_k_exp[[j]]-true_parameteres[[j]])^2))
}

#-------------------------------
# 5.3.2. LINES
#-------------------------------

graph_name = "Figure 5.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\textbf{{\\underline{\\hat{\\xi}}}}_{Whittle},\\textbf{{\\underline{\\xi}}})$")
xlab = "T"
main = ""
own_param = results_c_k_ar
r_param = results_c_k_exp
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param))
y_2 = colMeans((r_param))

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param)), 2, sd)
sd_r <- apply(((r_param)), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.4. MSE COEFFICIENTS FOR AR COEFFICIENT
################################################################################

# LOAD THE DATA
true_ar = readRDS(file = paste0(path,"true_ar.RData"))
fit_ar_ar = readRDS(file = paste0(path,"fit_ar_ar.RData"))
fit_ar_exp = readRDS(file = paste0(path,"fit_ar_exp.RData"))

#-------------------------------
# 5.4.2. LINES
#-------------------------------

graph_name = "Figure 6.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("MSE($\\hat{\\phi}_{1,Whittle},\\phi_1)$")
xlab = "T"
main = ""
true_param = true_ar
own_param = fit_ar_ar
r_param = fit_ar_exp
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param - true_param)^2)
y_2 = colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param - true_param)^2), 2, sd)
sd_r <- apply(((r_param - true_param)^2), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.5. MSE COEFFICIENTS FOR LAMBDA
################################################################################

# LOAD THE DATA
ar_lambda = readRDS(file = paste0(path,"ar_lambda.RData"))
exp_lambda = readRDS(file = paste0(path,"exp_lambda.RData"))

#-------------------------------
# 5.5.1. BOXPLOT
#-------------------------------

graph_name = "Figure 7.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{lambda}}_{MLE}$")
xlab = "T"
main = ""
y_1 = ar_lambda
y_2 = exp_lambda
abline_value = 1
abline_col = "red"

# LIMITS
lim_inf = min(y_1,y_2)
lim_sup = max(y_1,y_2)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,3,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

#-------------------------------
# 5.5.2. LINES
#-------------------------------

graph_name = "Figure 8.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{lambda}}_{MLE}$")
xlab = "T"
main =
true_param = 1
own_param = ar_lambda
r_param = exp_lambda
abline_value = 1
abline_col = "black"

# MSE
y_1 = colMeans((own_param))
y_2 = colMeans((r_param))

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param)), 2, sd)
sd_r <- apply(((r_param)), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.6. FREQUENCY DOMAIN GOODNESS OF FIT (FREQUENCY DOMAIN)
################################################################################

# LOAD THE DATA
ar_exp_1 = readRDS(file = paste0(path,"ar_exp_1.RData"))
exp_exp_1 = readRDS(file = paste0(path,"exp_exp_1.RData"))

#-------------------------------
# 5.6.1. BOXPLOT
#-------------------------------

graph_name = "Figure 9.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$p.value$")
xlab = "T"
main = ""
y_1 = ar_exp_1
y_2 = exp_exp_1
abline_value = 0.05
abline_col = "red"

# LIMITS
lim_inf = min(p_val_own_exp,p_val_r_exp)
lim_sup = max(p_val_own_exp,p_val_r_exp)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,3,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 1.5)
abline(h = abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

#-------------------------------
# 5.5.2. LINES
#-------------------------------

graph_name = "Figure 10.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("Non-rejected $H_0$")
xlab = "T"
main = ""
own_param = ar_exp_1
r_param = exp_exp_1
abline_value = 0.95
abline_col = "black"

# MSE
y_1 = colMeans((own_param>=0.05))
y_2 = colMeans((r_param>=0.05))

sd_own <- apply((own_param>=0.05), 2, sd)
sd_r <- apply((r_param>=0.05), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
dev.off()
