# fexpsmt
We are currently developing an R package that allows for the simulation of FARIMA and FEXP processes using their spectral density. This package is part of my Master's thesis project at the University of Geneva, where the advisor is Professor Dr. Davide La Vecchia, and the co-adviser is PhD student Manon Felix.

In the reproducibility folder one can find all the Monte Carlo experimnets that proove the correctness of this package. Each subfolder has two to three subfolders. The first can be "single fixed" which shows an example in depth for a specific parameter the respective process.The second can be "multiple fixed" which shows summary tables of the experiments for a significant sample of parameters, whicih vary from being close to non stationarity to being close of being white noise. The third can be "sampled" which simulates the time series with sampled parameters from a uniform distribution wihtin the bounds in which the process is stationary.

For the AR, MA, ARMA proceses, the experiments compare the simulations of fexpmst and stats packages, while for the FARIMA process they compare simulations of fexpmst and fracdiff package.
