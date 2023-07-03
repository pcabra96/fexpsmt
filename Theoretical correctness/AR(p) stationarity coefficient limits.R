library(latex2exp)
library(polynom)

num = 1


a =list()
a[[1]] = "kljaehfkjhdaf"
a[[2]] = rep(0,100)

num = 1
pr = polynomial(c(1,num))

p = 1

pr = pr*polynomial(c(1,num))
pr

################################################################################
# AR(1)
################################################################################

# Option 1
pr = polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 2
pr = polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# x-axis in (-1,1)

1

################################################################################
# AR(2)
################################################################################

# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 3
pr = polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

################################################################################
# AR(3)
################################################################################

# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 3
pr = polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 4
pr = polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

################################################################################
# AR(4)
################################################################################

# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 3
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 4
pr = polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 5
pr = polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

################################################################################
# AR(5)
################################################################################

# Option 1
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 2
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 3
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 4
pr = polynomial(c(1,num))*polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 5
pr = polynomial(c(1,num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef

# Option 6
pr = polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))*polynomial(c(1,-num))
pol_coef = -coefficients(pr)[-1]
pol_coef


generate_polynomial_coefficients <- function(order, num) {
  pol_coefs <- list()

  for (i in 1:order) {
    for (j in 1:3) {
      pr <- polynomial(rep(c(1, -1)[j], i)) * polynomial(c(1, num))
      pol_coef <- -coefficients(pr)[-1]
      pol_coefs[[paste0("Option", j, "_AR", i)]] <- pol_coef
    }
  }

  pol_coefs
}

order <- 3
num <- 2
coefficients <- generate_polynomial_coefficients(order, num)

# Print coefficients
for (key in names(coefficients)) {
  cat(key, ":", coefficients[[key]], "\n")
}

order = 2
for (i in 1:order){
  polynomial(c(1, num))^1
}







