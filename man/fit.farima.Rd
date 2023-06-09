% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.farima.R
\name{fit.farima}
\alias{fit.farima}
\title{Fit the FARIMA(p,d,q) parameters using Whittle estimation}
\usage{
fit.farima(y, p = 0, q = 0, d = 0)
}
\arguments{
\item{y}{Time Series}

\item{max.p}{Integer representing the order of the autoregressive component of the process}

\item{max.q}{Integer representing the order of the moving average component of the process}

\item{max.d}{1 if the process has long memory and 0 if it does not}
}
\value{
a list containing the estimated parameters of a FARIMA(p,d,q) process
}
\description{
Fit the FARIMA(p,d,q) parameters using Whittle estimation
}
\examples{
fit.farima(y, max.p = 1)
fit.farima(y, max.p = 1, max.q = 1, max.d = 1)
}
