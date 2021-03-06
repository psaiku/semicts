% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semicts.R
\name{rsemicts}
\alias{rsemicts}
\title{Generates a random sample from a semi-continuous distribution. Tested for
truncated normal (truncnorm), log-normal (lnorm), and gamma distributions are
supported for the continuous part of the distribution.}
\usage{
rsemicts(n, pzero = 0.5, r.func = NA, cts.density = "truncnorm",
  cts.param = list(a = 0, b = Inf, mean = 0, sd = 1))
}
\arguments{
\item{n}{Number of random variables to generate}

\item{pzero}{Point mass at 0}

\item{cts.density}{Name of a continuous density with support on the positive real line. Supported values: truncnorm (default), lnorm, and gamma}

\item{cts.param}{A list containing the parameters for cts.density (default: list(a=0, b=Inf, mean = 0, sd=1) for mean, standard deviation of truncated normal). For
log-normal, it should be an array containing meanlog, and sdlog of the distribution. For gamma,
an array of shape, and rate values must be supplied.}
}
\value{
An array of semi-continuous random variables.
}
\description{
Generates a random sample from a semi-continuous distribution. Tested for
truncated normal (truncnorm), log-normal (lnorm), and gamma distributions are
supported for the continuous part of the distribution.
}
\examples{
rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=list(meanlog=0, sdlog=1))
rsemicts(1000, pzero=0.6, cts.density="gamma", cts.param = list(shape=1, rate=1))
rsemicts(100, pzero=0.6, cts.density="truncnorm", cts.param = list(a=0, b=Inf, mean = 0, sd=1))
}
