% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dmvnrm}
\alias{dmvnrm}
\title{Calculate the multivariate normal density of a vector}
\usage{
dmvnrm(x, mean, sigma, logd = FALSE)
}
\arguments{
\item{x}{numeric vector containing the random numbers}

\item{mean}{numeric vector}

\item{sigma}{numeric covariance matrix}

\item{logd}{boolean indicating whether the log-transformed density
should be returned}
}
\description{
Calculate the multivariate normal density of a vector
}
