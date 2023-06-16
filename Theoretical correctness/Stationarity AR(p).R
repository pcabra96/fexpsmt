library(polynom)

AR_ORDER = 50
pr = polynomial(c(1))
for (i in 1:AR_ORDER) {
  pr = pr*polynomial(c(1,ifelse(runif(1)>=0.5,1,-1)*runif(n=1,min = -1,max = 1)))
}
pr
pol_coef = -coefficients(pr)[-1]

y = sim.farima(ar = pol_coef, T = 2^12)

plot(y[1:2^11], type = "l")

acf(y, lag.max = length(y)-1)

ifelse(runif(1)>=0.5,1,-1)
