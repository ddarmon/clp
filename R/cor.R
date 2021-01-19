#' Density of the Sample Pearson Correlation for a Bivariate Gaussian
#'
#' The density of the sample Pearson correlation for a bivariate Gaussian.
#'
#' @param x vector of quantiles
#' @param rho the population correlation coefficient
#' @param n the sample size
#'
#' @return A vector containing the density function evaluated at x.
#'
#' @references  Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
dcorr <- function(x, rho, n){
  prefactor <- function(r, rho, n){
    (n - 2)*(1 - rho^2)^((n-1)/2)*(1 - r^2)^((n-4)/2)/pi
  }

  integrand <- function(w, r, rho, n){
    1/(cosh(w) - rho*r)^(n-1)
  }

  definite.integral <- function(r, rho, n){
    integrate(integrand, lower = 0, upper = Inf, r = r, rho = rho, n = n)$value
  }
  definite.integral <- Vectorize(definite.integral, vectorize.args = 'r')

  I <- definite.integral(x, rho, n)

  fn <- prefactor(x, rho, n)*I

  return(fn)
}

dcorr <- Vectorize(dcorr)

#' Distribution Function of the Sample Pearson Correlation for a Bivariate Gaussian
#'
#' The distribution function of the sample Pearson correlation for a bivariate Gaussian.
#'
#' @param q vector of quantiles
#' @param rho the population correlation coefficient
#' @param n the sample size
#' @param lower.tail compute the lower tail (default) or
#'             upper tail of the confidence distribution.
#'
#' @return A vector containing the distribution function evaluated at x.
#'
#' @references  Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
pcorr <- function(q, rho, n, lower.tail = TRUE){
  if (lower.tail){
    integrate(dcorr, lower = -1, upper = q, rho = rho, n = n)$value
  }else{
    integrate(dcorr, lower = q, upper = 1, rho = rho, n = n)$value
  }

}

Vectorize(pcorr, vectorize.args = 'rho')

#' Confidence Functions for Pearson's Correlation Coefficient
#'
#' Confidence functions for Pearson's correlation coefficient
#' for a bivariate Gaussian.
#'
#' @param x a numeric vector
#' @param y a numeric vector
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#' @param exact whether the exact sampling distribution of the
#'              sample correlation coefficient (TRUE) or
#'              Fisher's Z-transformation (FALSE) should be
#'              used in constructing the confidence functions.
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for Pearson's correlation coefficient for a bivariate Gaussian, as well as
#'         the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
#' @examples
#' data(fat)
#'
#' fat <- fat[1:50, ] # Smaller sub-sample, to show exact versus
#'                    # Fisher's Z-transformation.
#'
#' # Using the exact sampling distribution of R
#' cor.conf(x = fat$body.fat, y = fat$weight, exact = TRUE)
#'
#' # Using Fisher's Z-transformation (to match cor.test())
#' cor.conf(x = fat$body.fat, y = fat$weight, exact = FALSE)
#'
#' @export cor.conf
cor.conf <- function(x, y, plot = TRUE, conf.level = 0.95, exact = FALSE){
  n <- length(x)

  if (n != length(y)){
    stop("x and y must be the same length.")
  }

  if (n <= 3){
    stop("n must be >= 4.")
  }

  if (n > 100 && exact == TRUE){
    warning("

            WARNING: n > 100.

            This can cause the integral needed for the sampling
            distribution of the correlation coefficient to become
            ill-behaved.

            Using Fisher's transformation instead.")
    exact <- FALSE
  }

  R <- cor(x, y)

  if (exact){ # Use exact sampling distribution of R
    pconf <- function(rho) {
      if (rho <= -1){
        return(0)
      }else if (rho >= 1){
        return(1)
      }else{
        return(pcorr(R, rho, n, lower.tail = FALSE))
      }
    }

    pconf <- Vectorize(pconf)

    dconf <- function(rho, dx = 1e-5) {
      if (rho <= -1 | rho >= 1){
        return(0)
      }else{
        return((pconf(rho + dx) - pconf(rho - dx))/(2*dx))
      }
    }

    dconf <- Vectorize(dconf)

    cconf <- function(rho) abs(2*pconf(rho) - 1)

    qconf <- function(p){
      fun.root <- function(x) pconf(x) - p

      return(uniroot(fun.root, interval = c(-1, 1))$root)
    }

    qconf <- Vectorize(qconf)
  }else{ # Use Fisher's Z-transformation
    Rt <- atanh(R)

    se <- 1/sqrt(n-3)

    pconf <- function(rho) pnorm((atanh(rho) - Rt)*sqrt(n - 3))

    dconf <- function(rho) dnorm((atanh(rho) - Rt)*sqrt(n - 3))*sqrt(n - 3)/(1 - rho^2)

    cconf <- function(rho) abs(2*pconf(rho) - 1)

    qconf <- function(p) tanh(Rt + qnorm(p)*se)
  }

  pcurve <- function(rho) 1 - cconf(rho)

  scurve <- function(rho) -log2(pcurve(rho))

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    plot.dconf(out, xlab = 'cor', n = 201)
    plot.cconf(out, conf.level = conf.level, xlab = 'cor', n = 201)
  }

  return(out)
}
