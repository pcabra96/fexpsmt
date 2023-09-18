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
SUBPROCESS = "fixed multiple"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/6. ",PROCESS,"/6.2. ",SUBPROCESS,"/")

# GENERAL
tite_size = 2.5
subtitle_size = 2
lab_size = 1.8
axis_text = 1.5
true_lab = 1.5

legend_size = 2
legend_line = 1
legend_line_style = 1.5

# LINES c(bottom, left, top, right))mar=c(5,5,6,2)
bottom_line = 5
left_line = 6
top_line = 6
right_line = 2

ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
T = 2^11

################################################################################
##----------------------------------------------------------------------------##
## 4. SIMULATION                                                              ##
##----------------------------------------------------------------------------##
################################################################################
dev.off()
for (j in 1:length(ar_coef_vec)) {
  ar_coef = ar_coef_vec[j]
  f_t_farima = farima.spectrum(ar = ar_coef,n.freq = T)

  mean_vec = matrix(0,nrow = 1000, ncol = 9)
  fitted_ar_vec = matrix(0,nrow = 1000, ncol = 9)
  num_ck = 1:9
  diff_vec = rep(0,length(num_ck))

  mhalfm <- (T-1) %/% 2L
  w <- 2*pi/T * (1:mhalfm)


  time_begin = Sys.time()
  for (i in num_ck) {
    coef_exp = fourier.series(log(f_t_farima), k = i)[["coef"]]
    f_t_fexp = fexp.spectrum(ck = coef_exp,n.freq = T)
    diff_vec[i] = trapz(w, abs(f_t_farima[1:(T/2-1)]-f_t_fexp[1:(T/2-1)]))
    for (sim in 1:1000) {
      y_exp = sim.fexp(ck = coef_exp, T = T)
      mean_vec[sim,i] = mean(y_exp) # OK
      fitted_ar_vec[sim,i] = fit.farima(y = y_exp,p = 1)[["ar"]] # OK
    }
  }
  time_end = Sys.time()
  time = time_end-time_begin
  print(time)

  # MEAN

  graph_name = paste0("Figure ",j,".pdf")
  pdf(file = paste0(path,graph_name), width = 8, height = 5)
  par(mfrow = c(1, 1), mar=c(bottom_line,3,1,right_line))
  plot(x = 1:9, y = diff_vec, type = "l", xlab = "p", ylab = "", cex.lab = lab_size*0.8, cex.axis = axis_text, ylim=c(0,0.4),xaxt = "n", col = "red")
  lines(x = 1:9, y = colMeans(abs(fitted_ar_vec-ar_coef)), type = "l", col = "blue")
  abline(h=0, col = "black")
  axis(1, at = 1:length(ar_coef_vec), labels = 1:length(ar_coef_vec))
  legend = c(TeX(paste0("$\\int_{0}^{\\pi}|f_{AR(\\phi_1=",ar_coef,")}(\\omega)-f_{EXP(p)}(\\omega)|\\partial \\omega$")),TeX("$\\bar{abs(\\hat{\\phi1_1}-\\phi_1)}$"),"0")
  legend("topright", legend = legend, col = c("red", "blue","black"), lty = 1,cex = legend_size*0.75,lwd=legend_line_style,seg.len=legend_line, bty="n")
  dev.off()
}




