% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop.R
\name{prop.conf}
\alias{prop.conf}
\title{Confidence Functions for One or Two Proportions}
\usage{
prop.conf(x, n, plot = TRUE, conf.level = 0.95)
}
\arguments{
\item{x}{the number of successes out of n trials (one proportion) or
a vector with the number of successes in two independent samples (two proportions).}

\item{n}{the number of trials (one proportion) or a vector with the number of trials
in two independent samples (two proportions)}

\item{plot}{whether to plot the confidence density and curve}

\item{conf.level}{the confidence level for the confidence interval indicated on the confidence curve}
}
\value{
A list containing the confidence functions pconf, dconf, cconf, and qconf
        for a single proportion or the difference of two proportions, as well as
        the P-curve and S-curve.
}
\description{
Confidence functions for a single binomial proportion based on the
mid P-value or the difference of two independent proportions based
on the score test statistic.
}
\examples{
# One proportion
 prop.conf(x = 0, n = 10)

# Two proportions
 prop.conf(x = c(1, 2), n = c(10, 20))

}
\references{
Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.

             Jan Hannig. "On generalized fiducial inference." Statistica Sinica (2009): 491-544.

             Olli Miettinen and Markku Nurminen. "Comparative analysis of two rates." Statistics in Medicine 4.2 (1985): 213-226.

             Markku Nurminen. "Confidence intervals for the difference and ratio of two binomial proportions." Biometrics 42 (1986): 675-676.
}
