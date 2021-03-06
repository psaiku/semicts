% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semicts.R
\name{summary.semicts}
\alias{summary.semicts}
\title{Prints the proportion of zeroes, and summary of the positive data in the
semicts object supplied.}
\usage{
\method{summary}{semicts}(obj, ...)
}
\arguments{
\item{obj}{A semicts object (for ex. returned from the rsemicts function)}
}
\description{
Prints the proportion of zeroes, and summary of the positive data in the
semicts object supplied.
}
\examples{
x <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
summary(x)
}
