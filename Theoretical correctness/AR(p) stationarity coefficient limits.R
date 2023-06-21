library(latex2exp)

num = 1

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
