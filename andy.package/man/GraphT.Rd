% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing_methods.R
\name{GraphT}
\alias{GraphT}
\title{Graph-based Two Sample Test}
\usage{
GraphT(x, ix, size, custom)
}
\arguments{
\item{x}{a matrix(n*d) of data points from the first and second group.
each row represent a data point;}

\item{ix}{a permutation of row indices of x.}

\item{size}{a vector of sample sizes;}

\item{custom}{the parameter which is customized. For this test,
custom means a vector used for the number of distance
and the threshold Q in graph-based test}
}
\value{
Statistics obtained after calculation
\item{Tr}{The graph-based two sample test statistics}
}
\description{
\code{GraphT} A graph-based two sample test for equal distributions
}
