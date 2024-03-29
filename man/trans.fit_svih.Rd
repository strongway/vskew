% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svi_calibration.R
\name{trans.fit_svih}
\alias{trans.fit_svih}
\title{fit SVI function with a known H and a known rho
Here we assume iv^H = SVI_Fun(Param)
The function will fit the 5 parameters, including exponent H, for the SVI function}
\usage{
trans.fit_svih(k, iv, H, rho, weights = 1)
}
\arguments{
\item{k}{log moneyness}

\item{iv}{iv}

\item{H}{H exponent}

\item{rho}{rho parameter in svi function}

\item{weights}{weights for individual rows (weighted optimization)}
}
\value{
svi 6 parameters
}
\description{
fit SVI function with a known H and a known rho
Here we assume iv^H = SVI_Fun(Param)
The function will fit the 5 parameters, including exponent H, for the SVI function
}
