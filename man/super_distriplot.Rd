% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScriptHomemadeFunctions.R
\name{super_distriplot}
\alias{super_distriplot}
\title{Creation of a multi-plot of variables distribution}
\usage{
super_distriplot(MYVARIABLES, GROUPS, breaks)
}
\arguments{
\item{MYVARIABLES}{A numeric vector or an array (up to 5) of numeric vectors (typically, variables as columns of a dataframe or matrix).}

\item{GROUPS}{A factor whose length matches the length of the variables of interest. The maximal number of groups/treatments to divide the observations
is limited to 4 (see 'note').}

\item{breaks}{A value, a vector of values, an algorithm etc. to set the width of the bars in the histograms (for more details on the possible input
values, please refer to the help page of \code{\link[graphics:hist]{hist}}).}
}
\value{
A vinwing window with as many histograms as there are combinations between the input variables and groups (up to 20).
}
\description{
The \code{super_distriplot()} function plots, in a single viewing window, the distribution of \code{MYVARIABLES} for all the \code{GROUPS} to
which they belong. It is very useful to quickly visualize if the \emph{among groups} \strong{normality assumption} is respected or not. The function
additionally plots the curve of the \emph{probability density function} of each sampled population, and that of a \strong{Normal} distribution.
}
\note{
For graphical convenience, the function cannot plot too many distributions in a single window. For that reason, only a maximum of 4 different
groups and 5 different variables is allowed as inputs in the function (giving a window with 20 plots). \cr
If ou wish to plot more variables and/or groups, please divide your data and run the function several times on each subset. \cr
Also, if you wish to plot a single variable divided in different groups, you need to specify \code{MYVARIABLES} that your \emph{single vector} (i.e.
your variable) is a dataframe. Otherwise, \code{super_distriplot()} will not be able to write a proper title to the histograms. If you want a proper
title, then use a synthax as follows: \code{... MYVARIABLES = as.data.frame(myvector) ...} or, for a \emph{single column}, \code{... MYVARIABLES =
as.data.frame(mydata$myvariable) ...} (see also 'examples').
}
\examples{
library(vegan)
data("dune")
data("dune.env")
super_distriplot(MYVARIABLES = dune[,3:6], GROUPS = dune.env[,4], breaks = 5)
# For a single variable (e.g. Bellis perennis):
super_distriplot(MYVARIABLES = as.data.frame(dune$Bellpere), GROUPS = dune.env[,4], breaks = 5)
}
