library(polynom)

AR_ORDER = 2
pr = polynomial(c(1))
for (i in 1:AR_ORDER) {
  pr = pr*polynomial(c(1,runif(n=1,min = 0.4,max = 0.8)))
}
pr
pol_coef = -coefficients(pr)[-1]
pol_coef

f_t = farima.spectrum(ar = pol_coef,n.freq = 2^11)
plot(x = 1:2^11, y = f_t, type = "l")
coef = fourier.series(log(f_t),k = 2)[["coef"]]
lines(x = 1:2^11, y = exp(fourier.series(log(f_t),k = 6)[["tfs"]]), type = "l", col = "red")

y = sim.farima(ar = pol_coef, T = 2^12)

plot(y[1:2^11], type = "l")

par(mfrow=c(1,1))
acf(y)
