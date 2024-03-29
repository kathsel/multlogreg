% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rtuvn}
\alias{rtuvn}
\title{Univariate Truncated Normal Distribution}
\usage{
rtuvn(n, mean, sd, lower, upper)
}
\arguments{
\item{n}{number of samples}

\item{mean}{mean}

\item{sd}{standard deviation}

\item{lower}{lower bound}

\item{upper}{upper bound}
}
\value{
vector of n samples
}
\description{
Generates n samples from the truncated normal distribution
 TN(mean, sd; lower, upper).
}
\examples{
rtuvn(n = 1, mean = 10, sd = 20, lower = 10, upper = 20)
rtuvn(n = 1, mean = 10, sd = 20, lower = 10, upper = Inf)

}
\references{
Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
    truncated multivariate normal and student-t distributions subject to
    linear inequality constraints. Journal of Statistical Theory and
    Practice, 9(4), 712-732.
}
