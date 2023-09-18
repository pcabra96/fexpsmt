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

################################################################################
# RANDOM POLYNOMIALS COEFFICIENTS
################################################################################

big_order = 50
num = 0.5
results = replicate(big_order, matrix(0,nrow = big_order+1, ncol = big_order), simplify = FALSE)
for (order in 2:big_order) {
  mat = matrix(0,nrow = big_order+1, ncol = big_order)
  for (j in 1:(order+1)) {
    a = runif(order,min = -num, max = num)
    for (po in 1:length(a)) {
      if (po==1) {
        pr=1
      }
      pr = pr*polynomial(c(1,a[po]))
      pol_coef = -coefficients(pr)[-1]
    }
    mat[j,1:order] = pol_coef
  }
  results[[order]] = mat
  names(results)[order] <- paste0("p = ",order)
}
results = results[-1]

################################################################################
# RESULTS
################################################################################

ar_order = 35
aux_matrix = results[[ar_order-1]]
n = 2^12
mhalfm <- (n-1) %/% 2L
w <- 2*pi/n * (1:mhalfm)
NCOL = 100

diff_matrix = matrix(0,nrow = (ar_order+1), ncol = NCOL)
for (row in 1:(ar_order+1)) {
  coef_ar = aux_matrix[row,1:ar_order]
  f_t_1 = farima.spectrum(ar = coef_ar, n.freq = n)
  coef_c_k = fourier.series(log(f_t_1), k = 100)[["coef"]]

  diff = rep(0,NCOL)
  for (i in 1:NCOL) {
    f_t_2 = fexp.spectrum(coef_c_k[1:i], n.freq = n)
    aux_diff = abs(f_t_1[1:length(w)]-f_t_2[1:length(w)])
    diff[i] = trapz(w, aux_diff)
  }
  diff_matrix[row,] = diff
}

boxplot(diff_matrix, ylim=c(0,5))
coL_means = colMeans(diff_matrix)
plot(coL_means, type = "l")



################################################################################
diff_matrix_2 = diff_matrix
diff_matrix = diff_matrix_2
diff_matrix = ifelse(diff_matrix<0.01,"0-0.01",
                     ifelse(diff_matrix<0.1 & diff_matrix>=0.01,"0.01-0.1",
                            ifelse(diff_matrix<1 & diff_matrix>=0.1,"0.1-1",
                                   ifelse(diff_matrix>=1 & diff_matrix<10,"1-10",
                                          ifelse(diff_matrix>=10 & diff_matrix<100,"10-100",diff_matrix)))))


# Create a custom color palette
col_palette_good = brewer.pal(7,"Blues")
col_palette_bad = brewer.pal(7,"Reds")
my_palette <- c("0-0.01" = col_palette_good[1], "0.01-0.1" = col_palette_good[2], "0.1-1" = col_palette_good[3],
                "1-10" = col_palette_bad[4], "10-100" = col_palette_bad[5],"100-1000"=col_palette_bad[6],"1000-Inf"=col_palette_bad[7])

# Convert diff_matrix to a data frame
df <- melt(diff_matrix)

# Assign colors to each type of string
df$colors <- factor(ifelse(df$value < 0.01, "0-0.01",
                           ifelse(df$value < 0.1 & df$value >= 0.01, "0.01-0.1",
                                  ifelse(df$value < 1 & df$value >= 0.1, "0.1-1",
                                         ifelse(df$value >= 1 & df$value < 10, "1-10",
                                                ifelse(df$value >= 10 & df$value < 100, "10-100",
                                                       ifelse(df$value >= 100 & df$value < 1000, "100-1000",
                                                              ifelse(df$value >= 1000, "1000-Inf", NA))))))),
                    levels = c("0-0.01", "0.01-0.1", "0.1-1", "1-10", "10-100","100-1000","1000-Inf"))

# Create a heatmap using ggplot2 with custom colors
ggplot(df, aes(x = Var2, y = Var1, fill = colors)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = my_palette, na.value = "gray") +
  labs(x = TeX(paste0("Order p of EXP(p)")), y = paste0("AR(",ar_order,")")) +
  ggtitle(TeX(paste0("$\\int_{0}^{\\pi}|f_{AR(",ar_order,")}(\\omega)-f_{EXP(p)}(\\omega)|d\\omega$"))) +
  theme_minimal()

# CHECK A SPECIFIC EXAPMLE
coef_ar = aux_matrix[8,1:ar_order]
f_t_1 = farima.spectrum(ar = coef_ar, n.freq = n)
coef_c_k = fourier.series(log(f_t_1), k = 100)[["coef"]]
diff = rep(0,50)
for (i in 1:50) {
  f_t_2 = fexp.spectrum(coef_c_k[1:i], n.freq = n)
  aux_diff = abs(f_t_1[1:length(w)]-f_t_2[1:length(w)])
  diff[i] = trapz(w, aux_diff)
}
par(mfrow=c(2,1))
plot(coef_ar[1:10], type = "o", main = "ar")
plot(coef_c_k[1:10], type = "o", main = "c_k")

par(mfrow=c(1,1))
plot(diff, type = "l", main = "Integral of ABS(AR(15) spectrum - EXP(p) expectrum)", xlab = "p", ylim=c(0,1))
abline(h=0.1, col ="red")







