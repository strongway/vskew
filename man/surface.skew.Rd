% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svi_calibration.R
\name{surface.skew}
\alias{surface.skew}
\title{plot skews from a surface}
\usage{
surface.skew(surf, dtes, strikes = seq(-0.2, 0.1, length.out = 300))
}
\arguments{
\item{surf}{surface of a given chain}

\item{dtes}{DTEs for skews}

\item{strikes}{range of strikes (or logmoneyness)}
}
\value{
estimated ivs (including skews, pdf etc. )
}
\description{
plot skews from a surface
}
