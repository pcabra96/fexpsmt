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
# 5.1. TIME (Figures: 1)
# 5.2. TIME SERIES AVERAGE (Figures: 2)
# 5.3. MSE FOR AR COEFFICIENT (Figures: 3)
# 5.4. MSE FOR MA COEFFICIENT (Figures: 4)
# 5.4. MSE FOR SHORT MEMORY PARAMETERS (Figures: 5)
# 5.5. MSE FOR LAMBDA COEFFICIENT (Figures: 6)
# 5.6. FREQUENCY DOMAIN GOODNESS OF FIT (Figures: 7)

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

# DO NOT TOUCH
active_path = dirname(rstudioapi::getActiveDocumentContext()$path)
path = paste0(active_path,"/")
ar_coef_vec = readRDS(file = paste0(path,"ar_coef_vec.RData"))
ma_coef_vec = readRDS(file = paste0(path,"ma_coef_vec.RData"))
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14
N_SIMULATIONS = 1000

################################################################################
##----------------------------------------------------------------------------##
## 4. VISUALIZATION PARAMETERS                                                ##
##----------------------------------------------------------------------------##
################################################################################

# GENERAL SIZES
tite_size = 2.5
subtitle_size = 2
lab_size = 1.5
axis_text = 1.5
true_lab = 1.5

# LEGEND DETAILS
graph_1 = "fepxmst"
graph_2 = "stats"
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
##----------------------------------------------------------------------------##
## 5. VISUALIZATION                                                           ##
##----------------------------------------------------------------------------##
################################################################################

################################################################################
# 5.1. TIME (Figures 1)
################################################################################

# LOAD THE DATA
time_own = readRDS(file = paste0(path,"time_own.RData"))
time_r = readRDS(file = paste0(path,"time_R.RData"))

#-------------------------------
# 5.1.1. LINES
#-------------------------------

graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = ""
y_1 = colMeans(time_own)
y_2 = colMeans(time_r)
abline_value = 0
abline_col = ""

# LIMITS
lim_inf = min(y_1,y_2)
lim_sup = max(y_1,y_2)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup+0.2), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.2. TIME SERIES AVERAGE (Figures: 2)
################################################################################

# LOAD THE DATA
average_own <- readRDS(file = paste0(path, "average_own.RData"))
average_r <- readRDS(file = paste0(path, "average_r.RData"))

#-------------------------------
# 5.2.2. LINES
#-------------------------------

graph_name = "Figure 2.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{mu}}$")
xlab = "T"
main = ""#paste0("MSE of $\\hat{mu}$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
true_param = 0
own_param = average_own
r_param = average_r
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans(own_param) #colMeans((own_param - true_param)^2)
y_2 = colMeans(r_param) #colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param )^2), 2, sd) # apply(((own_param - true_param)^2), 2, sd)
sd_r <- apply(((r_param)^2), 2, sd) # apply(((r_param - true_param)^2), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col, xpd = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

dev.off()

################################################################################
# 5.3. MSE FOR AR COEFFICIENT (Figures: 3)
################################################################################

# LOAD THE DATA
fit_own_phi <- readRDS(file = paste0(path, "phi_1_own.RData"))
fit_r_phi <- readRDS(file = paste0(path, "phi_1_r.RData"))

#-------------------------------
# 5.3.2. LINES
#-------------------------------

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\hat{phi}_{1,Whittle},\\phi_1)$")
xlab = "T"
main = ""
true_param = ar_coef_vec
own_param = fit_own_phi
r_param = fit_r_phi
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
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col, xpd = FALSE)
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
# 5.4. MSE FOR MA COEFFICIENT (Figures: 4)
################################################################################

# LOAD THE DATA
fit_own_theta = readRDS(file = paste0(path,"theta_1_own.RData"))
fit_r_theta = readRDS(file = paste0(path,"theta_1_r.RData"))

#-------------------------------
# 5.4.2. LINES
#-------------------------------

graph_name = "Figure 4.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\hat{theta}_{1,Whittle},\\theta_1)$")
xlab = "T"
main = ""
true_param = ma_coef_vec
own_param = fit_own_theta
r_param = fit_r_theta
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
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col, xpd = FALSE)
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
graph_name = "Figure 8.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.5. MSE FOR XI, SHORT PARAMETERS (Figures: 5)
################################################################################

#-------------------------------
# 5.4.2. LINES
#-------------------------------

# SAVE THE GRAPH
graph_name = "Figure 5.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\textbf{{\\underline{\\hat{\\xi}}}}_{Whittle},\\textbf{{\\underline{\\xi}}})$")
xlab = "T"
main = ""
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((fit_own_phi-ar_coef_vec)^2+(fit_own_theta-ma_coef_vec)^2)
y_2 = colMeans((fit_r_phi-ar_coef_vec)^2+(fit_r_theta-ma_coef_vec)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply((fit_own_phi-ar_coef_vec)^2+(fit_own_theta-ma_coef_vec)^2, 2, sd)
sd_r <- apply((fit_r_phi-ar_coef_vec)^2+(fit_r_theta-ma_coef_vec)^2, 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col, xpd = FALSE)
mtext(TeX(main), side = 3, line = -4, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

dev.off()


################################################################################
# 5.6. MSE FOR LAMBDA COEFFICIENT (Figures: 6)
################################################################################

# LOAD THE DATA
fit_own_exp = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp = readRDS(file = paste0(path,"lambda_r.RData"))

#-------------------------------
# 5.6.1. LINES
#-------------------------------

graph_name = "Figure 6.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\bar{\\hat{lambda}}_{MLE}$")
xlab = "T"
main = ""#paste0("Fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1}$")
true_param = 1
own_param = fit_own_exp
r_param = fit_r_exp
abline_value = 1
abline_col = "black"

# MSE
y_1 = colMeans((own_param))#colMeans((own_param - true_param)^2)
y_2 = colMeans((r_param))#colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply(((own_param )), 2, sd)
sd_r <- apply(((r_param)), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col, xpd = FALSE)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

dev.off()

################################################################################
# 5.7. FREQUENCY DOMAIN GOODNESS OF FIT (Figure 7)
################################################################################

# LOAD THE DATA
p_val_own_exp = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp = readRDS(file = paste0(path,"p.val_r.RData"))

# NAME
graph_name = "Figure 7.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("Non-rejected $H_0$ (%)")
xlab = "T"
main = ""#paste0("Percentage of non-rejected $H_0: \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
own_param = p_val_own_exp
r_param = p_val_r_exp
abline_value = 0.95
abline_col = "black"

# MSE
y_1 = colMeans((own_param>=0.05))
y_2 = colMeans((r_param>=0.05))

# STANDARD DEVIATION OF SQUARED ERROR
sd_own <- apply((own_param>=0.05), 2, sd)
sd_r <- apply((r_param>=0.05), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own, y_2 - sd_r)
lim_sup = max(y_1 + sd_own, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line, left_line, top_line, right_line), xpd = TRUE)

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend("topright", inset = c(inset_num, 0), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
mtext(TeX(main), side = 3, line = -1, outer = TRUE, cex=tite_size, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own[i], x1 = POWER[i], y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own[i], col = "blue")
  segments(x0 = POWER[i], y0 = y_2[i] - sd_r[i], x1 = POWER[i], y1 = y_2[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_r[i], col = "red")
}

dev.off()
