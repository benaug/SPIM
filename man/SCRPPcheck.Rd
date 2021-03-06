% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCRPPcheck.R
\name{SCRPPcheck}
\alias{SCRPPcheck}
\title{Calculate goodness of fit statistics for a model fit by mcmc.SCR. Currently won't calculate
#necessary statistics for individual-level p-values.}
\usage{
SCRPPcheck(data, posterior, use)
}
\arguments{
\item{data}{a list produced by simSCR2DNA or in the same format}

\item{posterior}{a posterior produced by mcmc.SCR}

\item{use}{a vector if iteration values at which the statistics should be calculated}
}
\description{
Calculate goodness of fit statistics for a model fit by mcmc.SCR
}
\author{
Ben Augustine
}
