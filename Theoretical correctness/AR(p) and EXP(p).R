library(latex2exp)

num = 0.7

################################################################################
# AR(1)
################################################################################
T = 2^12
# Option 1
pr = polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 2
pr = polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

################################################################################
# AR(2)
################################################################################

num = 0.8
T = 2^12
# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main), lag.max = 2^10)
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 3
pr = polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

################################################################################
# AR(3)
################################################################################

num = 0.9
# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 3
pr = polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 4
pr = polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

################################################################################
# AR(4)
################################################################################

# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 3
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 4
pr = polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")

# Option 5
pr = polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef
y = arima.sim(model = list(ar = pol_coef), n = T)
main = "$"
for (i in 1:length(pol_coef)) {
  main = paste0(main,"\\phi_",i," = ",pol_coef[i]," \\ ")
}
main = paste0(main,"$")
plot(y,type = "l", main = TeX(main))
acf(y, main = TeX(main))
f_t = farima.spectrum(ar=pol_coef,n.freq = T)

plot(x = 1:T,f_t, type = "l", main = TeX(main))
fourier_series = fourier.series(f_t = log(f_t),k = 2^10)
lines(fexp.spectrum(ck = fourier_series$coef,n.freq = T), type = "l", col = "red")
plot(fourier_series$coef[1:10], type = "l", ylab = TeX("$A_k$"), xlab = "k", main = "Cepstral coefficients")
abline(h=0,col = "red")
