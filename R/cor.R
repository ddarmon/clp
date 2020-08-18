#' Confidence Functions for Pearson's Correlation Coefficient
#'
#' Confidence functions for Pearson's correlation coefficient
#' via Fisher's Z-transformation.
#'
#' @param "x, y" numeric vectors of data values. x and y must have the same length.
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for Pearson's correlation coefficient via Fisher's Z-transformation, as well as
#'         the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' data(fat)
#' cor.conf(x = fat$body.fat, y = fat$weight)
#'
#' @export cor.conf
cor.conf <- function(x, y, plot = TRUE, conf.level = 0.95){
  n <- length(x)

  if (n != length(y)){
    stop("x and y must be the same length.")
  }

  if (n <= 3){
    stop("n must be >= 4.")
  }

  R <- cor(x, y)

  Rt <- atanh(R)

  se <- 1/sqrt(n-3)

  pconf <- function(rho) pnorm((atanh(rho) - Rt)*sqrt(n - 3))

  dconf <- function(rho) dnorm((atanh(rho) - Rt)*sqrt(n - 3))*sqrt(n - 3)/(1 - rho^2)

  cconf <- function(rho) abs(2*pconf(rho) - 1)

  qconf <- function(p) tanh(Rt + qnorm(p)*se)

  pcurve <- function(rho) 1 - cconf(rho)

  scurve <- function(rho) -log2(pcurve(rho))

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    plot.dconf(out, xlab = 'cor')
    plot.cconf(out, conf.level = conf.level, xlab = 'cor')
  }

  return(out)
}
