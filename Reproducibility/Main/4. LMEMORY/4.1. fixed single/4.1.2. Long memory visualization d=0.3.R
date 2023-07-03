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
# 5.3. MSE FOR d COEFFICIENT (Figures: 5,6)
# 5.4. MSE FOR LAMBDA COEFFICIENT, FREQUENCY DOMAIN (Figures: 9,10)
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT, FREQUENCY DOMAIN (Figures: 11, 12)

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

PROCESS = "LMEMORY"
SUBPROCESS = "fixed single"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/4. ",PROCESS,"/4.1. ",SUBPROCESS,"/")
d_coef_vec = 0.3
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
symbol = "d"

if (length(d_coef_vec)==1) {
  process_string = paste0("FARIMA(0,",d_coef_vec,",0)_t \\ or \\ FEXP(0,",d_coef_vec,"))")
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
axis_text = 1.5
true_lab = 1.5

graph_1 = paste0("$fexpmst_{FARIMA}$")
graph_2 = paste0("$fexpmst_{FEXP}$")
graph_3 = paste0("$fracdiff_{FARIMA}$")
legend = c(graph_1, graph_2, graph_3)
legend_size = 2
legend_line = 1
legend_line_style = 1.5

width_graphs = 1100
height_graphs = 550

# BOXPLOT c(bottom, left, top, right))mar=c(5,5,6,2)
bottom_box = 5
left_box = 6
top_box = 7.5
right_box = 2

# LINES c(bottom, left, top, right))mar=c(5,5,6,2)
bottom_line = 5
left_line = 6
top_line = 4
right_line = 2

########################################################
# 5.1. TIME (Figures: 1,2)
########################################################

# LOAD THE DATA
own_times_farima = readRDS(paste0(path,"own_times_farima.RData"))
own_times_fexp = readRDS(paste0(path,"own_times_fexp.RData"))
r_times = readRDS(paste0(path,"r_times.RData"))

#-------------------------------
# 5.1.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = paste0("Simulation time for $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
y_1 = own_times_farima
y_2 = own_times_fexp
y_3 = r_times
abline_value = 0
abline_col = ""

# LIMITS
lim_inf = min(y_1,y_2,y_3)
lim_sup = max(y_1,y_2,y_3)

# DISPLAY
par(mfrow=c(1,3), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_1), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_2), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
boxplot(y_3, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_3), cex.main = subtitle_size, line = 1)
mtext(TeX(main), side = 3, line = -5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 1.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.1.2. LINES
#-------------------------------

# PARAMETERS
ylab = "time (s)"
xlab = "T"
main = paste0("Average running time for $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
y_1 = colMeans(own_times_farima)
y_2 = colMeans(own_times_fexp)
y_3 = colMeans(r_times)
abline_value = 0
abline_col = ""

# LIMITS
min_lim = min(y_1, y_2, y_3)
max_lim = max(y_1, y_2, y_3)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, type = "o", ylim=c(min_lim,max_lim), col = "blue", ylab = ylab, xlab = xlab, labels = FALSE, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, type = "o", col = "green")
lines(x = POWER, y = y_3, type = "o", col = "red")
legend("topleft", legend = c(TeX(graph_1),TeX(graph_2), TeX(graph_3)), col = c("blue", "green", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
axis(1, at=POWER, labels = names)
axis(2)
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

graph_name = "Figure 2.png"
dev.print(device = png, filename = paste0(path,graph_name), width = 1800, height = 1100, res=200)
dev.off()
################################################################################
# 5.2. AVERAGE (Figures: 3,4)
################################################################################

# LOAD THE DATA
own_average_farima <- readRDS(file = paste0(path, "own_average_farima.RData"))
own_average_fexp <- readRDS(file = paste0(path, "own_average_fexp.RData"))
r_average <- readRDS(file = paste0(path, "r_average.RData"))

t.test(own_average_farima[,1],r_average[,1])
t.test(own_average_farima[,2],r_average[,2])
t.test(own_average_farima[,3],r_average[,3])
t.test(own_average_farima[,4],r_average[,4])
t.test(own_average_farima[,5],r_average[,5])
t.test(own_average_farima[,6],r_average[,6])
t.test(own_average_farima[,7],r_average[,7])
t.test(own_average_farima[,8],r_average[,8])

#-------------------------------
# 5.2.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{mu}$")
xlab = "T"
main = paste0("Fitted $\\mu$ of $\\",N_SIMULATIONS," \\ \\y_{",process_string,"_t}\\}_{t=1}^{T}$")
y_1 = own_average_farima
y_2 = own_average_fexp
y_3 = r_average
abline_value = 0
abline_col = "red"

# LIMITS
lim_inf = min(y_1,y_2, y_3)
lim_sup = max(y_1,y_2, y_3)

# DISPLAY
par(mfrow=c(1,3), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_1), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\mu$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_2), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\mu$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_3, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_3), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\mu$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

mtext(TeX(main), side = 3, line = -5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 3.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.2.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{mu}$ MSE")
xlab = "T"
main = paste0("MSE of $\\hat{mu}$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
true_param = 0
own_param_farima = own_average_farima
own_param_fexp = own_average_fexp
r_param = r_average
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param_farima - true_param)^2)
y_2 = colMeans((own_param_fexp - true_param)^2)
y_3 = colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own_farima <- apply(((own_param_farima - true_param)^2), 2, sd)
sd_own_fexp <- apply(((own_param_fexp - true_param)^2), 2, sd)
sd_r <- apply(((r_param - true_param)^2), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own_farima,y_2 - sd_own_fexp, y_2 - sd_r)
lim_sup = max(y_1 + sd_own_farima,y_2 + sd_own_fexp, y_2 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, type = "o", ylim=c(lim_inf,lim_sup), col = "blue", ylab = ylab, xlab = xlab, labels = FALSE, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, type = "o", col = "green")
lines(x = POWER, y = y_3, type = "o", col = "red")
legend("topleft", legend = c(TeX(graph_1),TeX(graph_2), TeX(graph_3)), col = c("blue", "green", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
axis(1, at=POWER, labels = names)
axis(2)
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i], y1 = y_1[i] + sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own_farima[i], col = "blue")

  segments(x0 = POWER[i], y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i], y1 = y_2[i] + sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_own_fexp[i], col = "green")

  segments(x0 = POWER[i], y0 = y_3[i] - sd_r[i], x1 = POWER[i], y1 = y_3[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
graph_name = "Figure 4.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.3. MSE FOR d COEFFICIENT (Figures: 5,6)
################################################################################

# LOAD THE DATA
own_long_param_farima <- readRDS(file = paste0(path, "own_long_param_farima.RData"))
own_long_param_fexp <- readRDS(file = paste0(path, "own_long_param_fexp.RData"))
r_long_param <- readRDS(file = paste0(path, "r_long_param.RData"))

t.test(own_long_param_farima[,1],r_long_param[,1])
t.test(own_long_param_farima[,2],r_long_param[,2])
t.test(own_long_param_farima[,3],r_long_param[,3])
t.test(own_long_param_farima[,4],r_long_param[,4])
t.test(own_long_param_farima[,5],r_long_param[,5])
t.test(own_long_param_farima[,6],r_long_param[,6])
t.test(own_long_param_farima[,7],r_long_param[,7])
t.test(own_long_param_farima[,8],r_long_param[,8])

#-------------------------------
# 5.3.1. BOXPLOT
#-------------------------------

if (length(d_coef_vec)==1) {
  # PARAMETERS
  ylab = TeX("$\\hat{d}_{MLE}$")
  xlab = "T"
  main = paste0("Fitted $\\d$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
  y_1 = own_long_param_farima
  y_2 = own_long_param_fexp
  y_3 = r_long_param
  abline_value = d_coef_vec
  abline_col = "red"

  # LIMITS
  lim_inf = min(y_1,y_2, y_3)
  lim_sup = max(y_1,y_2, y_3)

  # DISPLAY
  par(mfrow=c(1,3), mar=c(bottom_box,left_box,top_box,right_box))

  # PLOT
  boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = TeX(graph_1), cex.main = subtitle_size, line = 1)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\d$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

  boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = TeX(graph_2), cex.main = subtitle_size, line = 1)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\d$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

  boxplot(y_3, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
  title(main = TeX(graph_3), cex.main = subtitle_size, line = 1)
  abline(h = abline_value, col = abline_col)
  legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\d$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

  mtext(TeX(main), side = 3, line = -5, outer = TRUE, cex=tite_size, font = 2)

  # SAVE THE GRAPH
  graph_name = "Figure 5.png"
  dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
  dev.off()
}

#-------------------------------
# 5.3.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{d}_{MLE}$ MSE")
xlab = "T"
main = paste0("MSE of $\\hat{d}$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
true_param = d_coef_vec
own_param_farima = own_long_param_farima
own_param_fexp = own_long_param_fexp
r_param = r_long_param
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param_farima - true_param)^2)
y_2 = colMeans((own_param_fexp - true_param)^2)
y_3 = colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own_farima <- apply(((own_param_farima - true_param)^2), 2, sd)
sd_own_fexp <- apply(((own_param_fexp - true_param)^2), 2, sd)
sd_r <- apply(((r_param - true_param)^2), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own_farima, y_2 - sd_own_fexp, y_3 - sd_r)
lim_sup = max(y_1 + sd_own_farima, y_2 + sd_own_fexp, y_3 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, type = "o", ylim=c(lim_inf,lim_sup), col = "blue", ylab = ylab, xlab = xlab, labels = FALSE, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, type = "o", col = "green")
lines(x = POWER, y = y_3, type = "o", col = "red")
legend("topright", legend = c(TeX(graph_1),TeX(graph_2), TeX(graph_3)), col = c("blue", "green", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
axis(1, at=POWER, labels = names)
axis(2)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i], y1 = y_1[i] + sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own_farima[i], col = "blue")

  segments(x0 = POWER[i], y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i], y1 = y_2[i] + sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_own_fexp[i], col = "green")

  segments(x0 = POWER[i], y0 = y_3[i] - sd_r[i], x1 = POWER[i], y1 = y_3[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
graph_name = "Figure 6.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.4. MSE FOR LAMBDA COEFFICIENT, FREQUENCY DOMAIN (Figures: 9,10)
################################################################################

# LOAD THE DATA
own_lambda_farima = readRDS(file = paste0(path,"own_lambda_farima.RData"))
own_lambda_fexp = readRDS(file = paste0(path,"own_lambda_fexp.RData"))
r_lambda = readRDS(file = paste0(path,"r_lambda.RData"))

t.test(own_lambda_farima[,1],r_lambda[,1])
t.test(own_lambda_farima[,2],r_lambda[,2])
t.test(own_lambda_farima[,3],r_lambda[,3])
t.test(own_lambda_farima[,4],r_lambda[,4])
t.test(own_lambda_farima[,5],r_lambda[,5])
t.test(own_lambda_farima[,6],r_lambda[,6])
t.test(own_lambda_farima[,7],r_lambda[,7])
t.test(own_lambda_farima[,8],r_lambda[,8])

#-------------------------------
# 5.4.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{lambda}_{MLE}$")
xlab = "T"
main = paste0("Fitted $\\lambda$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
y_1 = own_lambda_farima
y_2 = own_lambda_fexp
y_3 = r_lambda
abline_value = 1
abline_col = "red"

# LIMITS
lim_inf = min(y_1,y_2, y_3)
lim_sup = max(y_1,y_2, y_3)

# DISPLAY
par(mfrow=c(1,3), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_1), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_2), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_3, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_3), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("True $\\lambda$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

mtext(TeX(main), side = 3, line = -5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 7.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.4.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("$\\hat{lambda}_{MLE}$ MSE")
xlab = "T"
main = paste0("MSE of $\\hat{lambda}$ with $\\",N_SIMULATIONS," \\ \\{y_{",process_string,"_t}\\}_{t=1}^{T}$")
true_param = 1
own_param_farima = own_lambda_farima
own_param_fexp = own_lambda_fexp
r_param = r_lambda
abline_value = 0
abline_col = "black"

# MSE
y_1 = colMeans((own_param_farima - true_param)^2)
y_2 = colMeans((own_param_fexp - true_param)^2)
y_3 = colMeans((r_param - true_param)^2)

# STANDARD DEVIATION OF SQUARED ERROR
sd_own_farima <- apply(((own_param_farima - true_param)^2), 2, sd)
sd_own_fexp <- apply(((own_param_fexp - true_param)^2), 2, sd)
sd_r <- apply(((r_param - true_param)^2), 2, sd)

# LIMITS
lim_inf = min(y_1 - sd_own_farima, y_2 - sd_own_fexp, y_3 - sd_r)
lim_sup = max(y_1 + sd_own_farima, y_2 + sd_own_fexp, y_3 + sd_r)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, type = "o", ylim=c(lim_inf,lim_sup), col = "blue", ylab = ylab, xlab = xlab, labels = FALSE, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, type = "o", col = "green")
lines(x = POWER, y = y_3, type = "o", col = "red")
legend("topright", legend = c(TeX(graph_1),TeX(graph_2), TeX(graph_3)), col = c("blue", "green", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
axis(1, at=POWER, labels = names)
axis(2)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

# PLOT STANDARD DEVIATION
for (i in 1:length(POWER)) {
  segments(x0 = POWER[i], y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i], y1 = y_1[i] + sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] - sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] - sd_own_farima[i], col = "blue")
  segments(x0 = POWER[i] - 0.1, y0 = y_1[i] + sd_own_farima[i], x1 = POWER[i] + 0.1, y1 = y_1[i] + sd_own_farima[i], col = "blue")

  segments(x0 = POWER[i], y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i], y1 = y_2[i] + sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] - sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] - sd_own_fexp[i], col = "green")
  segments(x0 = POWER[i] - 0.1, y0 = y_2[i] + sd_own_fexp[i], x1 = POWER[i] + 0.1, y1 = y_2[i] + sd_own_fexp[i], col = "green")

  segments(x0 = POWER[i], y0 = y_3[i] - sd_r[i], x1 = POWER[i], y1 = y_3[i] + sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] - sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] - sd_r[i], col = "red")
  segments(x0 = POWER[i] - 0.1, y0 = y_3[i] + sd_r[i], x1 = POWER[i] + 0.1, y1 = y_3[i] + sd_r[i], col = "red")
}

# SAVE THE GRAPH
graph_name = "Figure 8.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT, FREQUENCY DOMAIN (Figures: 11, 12)
################################################################################

# LOAD THE DATA
own_exp_1_farima = readRDS(file = paste0(path,"own_exp_1_farima.RData"))
own_exp_1_fexp = readRDS(file = paste0(path,"own_exp_1_fexp.RData"))
r_exp_1 = readRDS(file = paste0(path,"r_exp_1.RData"))

t.test(own_exp_1_farima[,1],r_exp_1[,1])
t.test(own_exp_1_farima[,2],r_exp_1[,2])
t.test(own_exp_1_farima[,3],r_exp_1[,3])
t.test(own_exp_1_farima[,4],r_exp_1[,4])
t.test(own_exp_1_farima[,5],r_exp_1[,5])
t.test(own_exp_1_farima[,6],r_exp_1[,6])
t.test(own_exp_1_farima[,7],r_exp_1[,7])
t.test(own_exp_1_farima[,8],r_exp_1[,8])

#-------------------------------
# 5.5.1. BOXPLOT
#-------------------------------

# PARAMETERS
ylab = TeX("$p.value$")
xlab = "T"
main = paste0("$H_0: \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
y_1 = own_exp_1_farima
y_2 = own_exp_1_fexp
y_3 = r_exp_1
abline_value = 0.05
abline_col = "red"

# LIMITS
lim_inf = min(y_1,y_2, y_3)
lim_sup = max(y_1,y_2, y_3)

# DISPLAY
par(mfrow=c(1,3), mar=c(bottom_box,left_box,top_box,right_box))

# PLOT
boxplot(y_1, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_1), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("$\\alpha=0.05$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_2, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_2), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("$\\alpha=0.05$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

boxplot(y_3, ylim=c(lim_inf,lim_sup), names = names, ylab = ylab, xlab = xlab, cex.lab = lab_size, cex.axis = axis_text)
title(main = TeX(graph_3), cex.main = subtitle_size, line = 1)
abline(h = abline_value, col = abline_col)
legend(3,(lim_inf+(lim_sup-lim_inf)/6), legend= TeX("$\\alpha=0.05$"),col = "red", lty=1:2, cex=true_lab, bty="n",lwd=legend_line_style,seg.len=legend_line)

mtext(TeX(main), side = 3, line = -5, outer = TRUE, cex=tite_size, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 9.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()

#-------------------------------
# 5.5.2. LINES
#-------------------------------

# PARAMETERS
ylab = TeX("Non-rejected $H_0$")
xlab = "T"
main = paste0("Percentage of non-rejected $H_0: \\ \\{I^*_{",process_string,"}(\\omega_k)\\}_{k=1}^{T-1} \\sim exp(\\lambda=1)$")
own_param_farima = own_exp_1_farima
own_param_fexp = own_exp_1_fexp
r_param = r_exp_1
abline_value = 0.95
abline_col = "black"

# MSE
y_1 = colMeans((own_param_farima>=0.05))
y_2 = colMeans((own_param_fexp>=0.05))
y_3 = colMeans((r_param>=0.05))

# STANDARD DEVIATION OF SQUARED ERROR
sd_own_farima <- apply(((own_param_farima>=0.05)), 2, sd)
sd_own_fexp <- apply(((own_param_fexp>=0.05)), 2, sd)
sd_r <- apply(((r_param>=0.05)), 2, sd)

# LIMITS
lim_inf = min(0.9)
lim_sup = max(1)

# DISPLAY
par(mfrow = c(1, 1), mar=c(bottom_line,left_line,top_line,right_line))

# PLOT
plot(x = POWER, y = y_1, type = "o", ylim=c(lim_inf,lim_sup), col = "blue", ylab = ylab, xlab = xlab, labels = FALSE, cex.lab = lab_size, cex.axis = axis_text)
lines(x = POWER, y = y_2, type = "o", col = "green")
lines(x = POWER, y = y_3, type = "o", col = "red")
legend("topright", legend = c(TeX(graph_1),TeX(graph_2), TeX(graph_3)), col = c("blue", "green", "red"), lty = 1, cex = legend_size,lwd=legend_line_style,seg.len=legend_line, bty="n")
axis(1, at=POWER, labels = names)
axis(2)
abline(h=abline_value, col = abline_col)
mtext(TeX(main), side = 3, line = -2.5, outer = TRUE,cex=1.5, font = 2)

# SAVE THE GRAPH
graph_name = "Figure 10.png"
dev.print(device = png, filename = paste0(path,graph_name), width = width_graphs, height = height_graphs)
dev.off()
