# SIMULATION PARAMETERS
NUM_SIM = 1000
POWER = 7:14
fitted_results = matrix(0,nrow = NUM_SIM, ncol = length(POWER))
c_k = 1

# SIMULATION
for (j in 1:length(7:14)) {
  T = 2^(POWER[j])
  for (sim in 1:NUM_SIM) {
    y = sim.fexp(c_k, T = T)
    fitted_results[sim,j] = fit.fexp(y, p = 1)[["c_k"]]
  }
}

# RESULTS
sum((fitted_results-c_k)^2)/NUM_SIM
boxplot(fitted_results, names=c("128", "256", "512", "1024", "2048", "4096", "8192","16384"),
                        xlab = expression('Simulated {Y' [t]~'} length T'), ylab = expression('Fitted A'[1]))
abline(h=c_k, col = "red")

mse = apply(fitted_results, 2, function(x) mean((x - c_k)^2))
