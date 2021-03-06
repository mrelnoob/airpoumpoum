% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScriptHomemadeFunctions.R
\name{se.m}
\alias{se.m}
\title{Standard Error of a mean}
\usage{
se.m(x)
}
\arguments{
\item{x}{a vector}
}
\value{
The standard error of the input vector's mean.
}
\description{
Computes the standard error of a sample mean. Useful to add error bars on plots or to compute intervals of confidence.
}
\examples{
myvector <- rnorm(n = 100, mean = 48, sd = 16)
se.m(x = myvector)
}
