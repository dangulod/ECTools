% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{k5}
\alias{k5}
\title{The fifth central moment}
\usage{
k5(x)
}
\arguments{
\item{loss}{vector with input}

\item{scale}{scale parameter}
}
\value{
High-order moments are moments beyond 4th-order moments. As with variance, skewness, and kurtosis, these are higher-order statistics,
involving non-linear combinations of the data, and can be used for description or estimation of further shape parameters. The higher the moment,
the harder it is to estimate, in the sense that larger samples are required in order to obtain estimates of similar quality. This is due to the
excess degrees of freedom consumed by the higher orders. Further, they can be subtle to interpret, often being most easily understood in terms of
lower order moments – compare the higher derivatives of jerk and jounce in physics. For example, just as the 4th-order moment (kurtosis) can be
interpreted as "relative importance of tails versus shoulders in causing dispersion" (for a given dispersion, high kurtosis corresponds to heavy
tails, while low kurtosis corresponds to broad shoulders), the 5th-order moment can be interpreted as measuring "relative importance of tails
versus center (mode, shoulders) in causing skew" (for a given skew, high 5th moment corresponds to heavy tail and little movement of mode, while
low 5th moment corresponds to more change in shoulders).
}
\description{
The fifth central moment
}
\references{
https://en.wikipedia.org/wiki/Moment_(mathematics)
}
