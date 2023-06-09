% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectral.density.fexp.R
\name{spectral.density.fexp}
\alias{spectral.density.fexp}
\title{Spectral density of an arbitrary FEXP process.}
\usage{
spectral.density.fexp(eta = c(), p = 0, d_long = 0, x)
}
\arguments{
\item{eta}{Vector with the parameters}

\item{p}{Integer representing the order of the FEXP process}

\item{d_long}{Integer with zero or one representing if one considers the long memory component}

\item{x}{Frequencies to be used}
}
\value{
Spectral density of the arbitrary FARIMA process
}
\description{
It is used to fit the simulated time series.
}
\examples{
spectral.density.farima(eta = c(0.5,0.5,0.3), p = 2, d_long=1, x)
}
