% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop.R
\name{matchedprop.conf}
\alias{matchedprop.conf}
\title{Confidence Functions for the Difference of Two Matched Proportions}
\usage{
matchedprop.conf(b, c, n, plot = TRUE, conf.level = 0.95)
}
\arguments{
\item{b}{the number of trials where the first outcome is a failure
and the second outcome is a success}

\item{c}{the number of trials where the first outcome is a success
and the second outcome is a failure}

\item{n}{the number of trials}

\item{plot}{whether to plot the confidence density and curve}

\item{conf.level}{the confidence level for the confidence interval indicated on the confidence curve}
}
\value{
A list containing the confidence functions pconf, dconf, cconf, and qconf
        for the difference of two matched proportions, as well as the P-curve and S-curve.
}
\description{
Confidence functions for the difference of two matched proportions
based on the score statistic.
}
\examples{
matchedprop.conf(b = 5, c = 10, n = 100)

}
\references{
Toshiro Tango. "Equivalence test and confidence interval for the difference in proportions for the paired‐sample design." Statistics in Medicine 17.8 (1998): 891-908.
}
