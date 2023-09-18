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

# TO TOUCH
PROCESS = "FARIMA"
SUBPROCESS = "sampled"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/5. ",PROCESS,"/5.2. ",SUBPROCESS,"/")
ar_coef_vec = readRDS(file = paste0(path,"ar_coef_vec.RData"))
ma_coef_vec = readRDS(file = paste0(path,"ma_coef_vec.RData"))
d_coef_vec = readRDS(file = paste0(path,"d_coef_vec.RData"))

# DO NOT TOUCH
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
POWER = 7:14
N_SIMULATIONS = 1000

if (length(ar_coef_vec)==1) {
  process_string = paste0(PROCESS,"(phi_1=",ar_coef_vec," , d=",d_coef_vec,", theta_1=",ma_coef_vec,")")
}else{
  process_string = paste0(PROCESS,"(phi_1,theta_1 \\sim u(\\[-0.9,0.9\\]))")
}

################################################################################
##----------------------------------------------------------------------------##
## 4. VISUALIZATION PARAMETERS                                                ##
##----------------------------------------------------------------------------##
################################################################################

# GENERAL
tite_size = 2.5
subtitle_size = 2
lab_size = 1.8
axis_text = 1.2
true_lab = 1.2

graph_1 = "fepxmst"
graph_2 = "fracdiff"
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

# LOAD THE DATA
time_own = readRDS(file = paste0(path,"time_own.RData"))
time_r = readRDS(file = paste0(path,"time_R.RData"))

#-------------------------------
# 5.1.1. BOXPLOT
#-------------------------------

graph_name = "Figure 1.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = ""
y_1 = time_own
y_2 = time_r
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
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

#-------------------------------
# 5.1.2. LINES
#-------------------------------
dev.off()
graph_name = "Figure 2.pdf"
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
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.2. Time Series average
################################################################################

graph_name = "Figure 3.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# LOAD THE DATA
average_own <- readRDS(file = paste0(path, "average_own.RData"))
average_r <- readRDS(file = paste0(path, "average_r.RData"))

t.test(average_own[,1],average_r[,1])
t.test(average_own[,2],average_r[,2])
t.test(average_own[,3],average_r[,3])
t.test(average_own[,4],average_r[,4])
t.test(average_own[,5],average_r[,5])
t.test(average_own[,6],average_r[,6])
t.test(average_own[,7],average_r[,7])
t.test(average_own[,8],average_r[,8])

#-------------------------------
# 5.2.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{mu}$ MSE")
xlab = "T"
main = ""
y_1 = average_own
y_2 = average_r
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
own_param = average_own
r_param = average_r
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
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
# 5.3. MSE COEFFICIENTS FOR AR
################################################################################

# LOAD THE DATA
fit_own_phi <- readRDS(file = paste0(path, "phi_1_own.RData"))
fit_r_phi <- readRDS(file = paste0(path, "phi_1_r.RData"))

#-------------------------------
# 5.3.1. BOXPLOT
#-------------------------------

if (length(ar_coef_vec)==1) {
  print(1)
  # PARAMETERS
  ylab = TeX("$\\hat{phi}_{1,\\ Whittle}$")
  xlab = "T"
  main = paste0("Fitted $\\phi_1$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
  y_1 = fit_own_phi
  y_2 = fit_r_phi
  abline_value = ar_coef_vec
  abline_col = "red"

  # LIMITS
  lim_inf = min(y_1,y_2)
  lim_sup = max(y_1,y_2)

  # DISPLAY
  par(mfrow=c(1,2), mar=c(bottom_box,left_box,top_box,right_box))

  # PLOT
  boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_1, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\phi_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_2, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\phi_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

  # SAVE THE GRAPH
  graph_name = "Figure 5.png"
  dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
  dev.off()
}

#-------------------------------
# 5.3.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{phi}_{1,\\ Whittle}$ MSE")
xlab = "T"
main = paste0("MSE of $\\hat{phi}_1$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
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
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
graph_name = "Figure 6.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.4. MSE COEFFICIENTS FOR MA COEFFICIENT
################################################################################

# LOAD THE DATA
fit_own_theta = readRDS(file = paste0(path,"theta_1_own.RData"))
fit_r_theta = readRDS(file = paste0(path,"theta_1_r.RData"))

t.test(fit_own_theta[,1],fit_r_theta[,1])
t.test(fit_own_theta[,2],fit_r_theta[,2])
t.test(fit_own_theta[,3],fit_r_theta[,3])
t.test(fit_own_theta[,4],fit_r_theta[,4])
t.test(fit_own_theta[,5],fit_r_theta[,5])
t.test(fit_own_theta[,6],fit_r_theta[,6])
t.test(fit_own_theta[,7],fit_r_theta[,7])
t.test(fit_own_theta[,8],fit_r_theta[,8])

#-------------------------------
# 5.4.1. BOXPLOT
#-------------------------------
if (length(ma_coef_vec)==1) {
  # PARAMETERS
  ylab = TeX("$\\hat{theta}_{1,\\ Whittle}$")
  xlab = "T"
  main = paste0("Fitted $\\theta_1$ with $\\",N_SIMULATIONS," \\ \\{y_{ARMA(phi_1=",ar_coef_vec,",theta_1 = ",ma_coef_vec,")_t}\\}_{t=1}^{T}$")
  y_1 = fit_own_theta
  y_2 = fit_r_theta
  abline_value = ma_coef_vec
  abline_col = "red"

  # LIMITS
  lim_inf = min(y_1,y_2)
  lim_sup = max(y_1,y_2)

  # DISPLAY
  par(mfrow=c(1,2), mar=c(bottom_box,left_box,top_box,right_box))

  # PLOT
  boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_1, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\theta_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_2, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\theta_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

  # SAVE THE GRAPH
  graph_name = "Figure 7.png"
  dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
  dev.off()
}

#-------------------------------
# 5.4.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{theta}_{1,\\ Whittle}$ $MSE")
xlab = "T"
main = paste0("MSE of $\\hat{theta}_1$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
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
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
graph_name = "Figure 8.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.5. MSE COEFFICIENTS FOR d COEFFICIENT
################################################################################

# LOAD THE DATA
d_own = readRDS(file = paste0(path,"d_own.RData"))
d_r = readRDS(file = paste0(path,"d_r.RData"))

t.test(d_own[,1],d_r[,1])
t.test(d_own[,2],d_r[,2])
t.test(d_own[,3],d_r[,3])
t.test(d_own[,4],d_r[,4])
t.test(d_own[,5],d_r[,5])
t.test(d_own[,6],d_r[,6])
t.test(d_own[,7],d_r[,7])
t.test(d_own[,8],d_r[,8])

#-------------------------------
# 5.4.1. BOXPLOT
#-------------------------------
if (length(d_coef_vec)==1) {
  # PARAMETERS
  ylab = TeX("$\\hat{d}_{\\ Whittle}$")
  xlab = "T"
  main = paste0("Fitted $\\d$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
  y_1 = d_own
  y_2 = d_r
  abline_value = d_coef_vec
  abline_col = "red"

  # LIMITS
  lim_inf = min(y_1,y_2)
  lim_sup = max(y_1,y_2)

  # DISPLAY
  par(mfrow=c(1,2), mar=c(bottom_box,left_box,top_box,right_box))

  # PLOT
  boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_1, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\theta_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = graph_2, cex.main = subtitle_size, line = 0.5)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\theta_1$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
  mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

  # SAVE THE GRAPH
  graph_name = "Figure 9.png"
  dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
  dev.off()
}

#-------------------------------
# 5.4.2. LINES
#-------------------------------

graph_name = "Figure 10.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\hat{d}_{Whittle},d)$")
xlab = "T"
main = ""
true_param = d_coef_vec
own_param = d_own
r_param = d_r
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
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
fit_own_exp = readRDS(file = paste0(path,"lambda_own.RData"))
fit_r_exp = readRDS(file = paste0(path,"lambda_r.RData"))

t.test(fit_own_exp[,1],fit_r_exp[,1])
t.test(fit_own_exp[,2],fit_r_exp[,2])
t.test(fit_own_exp[,3],fit_r_exp[,3])
t.test(fit_own_exp[,4],fit_r_exp[,4])
t.test(fit_own_exp[,5],fit_r_exp[,5])
t.test(fit_own_exp[,6],fit_r_exp[,6])
t.test(fit_own_exp[,7],fit_r_exp[,7])
t.test(fit_own_exp[,8],fit_r_exp[,8])

#-------------------------------
# 5.5.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{lambda}_{MLE}$")
xlab = "T"
main = paste0("Fitted $\\lambda$ for $\\",N_SIMULATIONS," \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1}$")
y_1 = fit_own_exp
y_2 = fit_r_exp
abline_value = 1
abline_col = "red"

# LIMITS
lim_inf = min(fit_own_exp,fit_r_exp)
lim_sup = max(fit_own_exp,fit_r_exp)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 0.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 0.5)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 11.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.5.2. LINES
#-------------------------------

graph_name = "Figure 12.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$\\hat{lambda}_{MLE})$")
xlab = "T"
main = ""
true_param = 1
own_param = fit_own_exp
r_param = fit_r_exp
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
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
graph_name = "Figure 12.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.6. FREQUENCY DOMAIN GOODNESS OF FIT (FREQUENCY DOMAIN)
################################################################################

# LOAD THE DATA
p_val_own_exp = readRDS(file = paste0(path,"p.val_own.RData"))
p_val_r_exp = readRDS(file = paste0(path,"p.val_r.RData"))

t.test(p_val_own_exp[,1],p_val_r_exp[,1])
t.test(p_val_own_exp[,2],p_val_r_exp[,2])
t.test(p_val_own_exp[,3],p_val_r_exp[,3])
t.test(p_val_own_exp[,4],p_val_r_exp[,4])
t.test(p_val_own_exp[,5],p_val_r_exp[,5])
t.test(p_val_own_exp[,6],p_val_r_exp[,6])
t.test(p_val_own_exp[,7],p_val_r_exp[,7])
t.test(p_val_own_exp[,8],p_val_r_exp[,8])

#-------------------------------
# 5.6.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$p.value$")
xlab = "T"
main = paste0("$H_0: \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
y_1 = p_val_own_exp
y_2 = p_val_r_exp
abline_value = 0.05
abline_col = "red"

# LIMITS
lim_inf = min(p_val_own_exp,p_val_r_exp)
lim_sup = max(p_val_own_exp,p_val_r_exp)

# DISPLAY
par(mfrow=c(1,2), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_1, cex.main = subtitle_size, line = 0.5)
abline(h = abline_value, col = abline_col)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = graph_2, cex.main = subtitle_size, line = 0.5)
abline(h = abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -3.5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 13.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.5.2. LINES
#-------------------------------

graph_name = "Figure 14.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("Non-rejected $H_0$")
xlab = "T"
main = ""
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
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend(11,(lim_sup+0.01), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
abline(h=abline_value, col = abline_col)
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

# SAVE THE GRAPH
dev.off()

################################################################################
# 5.7. MSE COEFFICIENTS FOR XI
################################################################################

#-------------------------------
# 5.4.2. LINES
#-------------------------------

graph_name = "Figure 15.pdf"
pdf(file = paste0(path,graph_name), width = 8, height = 5) # The height of the plot in inches

# PARAMETERS
ylab = TeX("$MSE(\\textbf{{\\underline{\\hat{\\xi}}}}_{Whittle},\\textbf{{\\underline{\\xi}}})$")
xlab = "T"
main = ""
true_param = ma_coef_vec
own_param = fit_own_theta
r_param = fit_r_theta
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
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,1,right_line))

# PLOT
plot(x = POWER, y = y_1, col = "blue", type = "o", ylim=c(lim_inf,lim_sup), ylab = ylab, labels = FALSE, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, col = "red", type = "o")
axis(1, at=POWER, labels = names)
axis(2)
legend(11,(lim_sup), legend = legend, col = c("blue", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n",)
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
