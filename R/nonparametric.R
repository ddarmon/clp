#' Confidence Functions for the Median via the Sign Test
#'
#' Confidence functions for the median of a continuous distribution based on the
#' sign test.
#'
#' @param x a vector of observations
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the median, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Myles Hollander, Douglas A. Wolfe, and Eric Chicken. Nonparametric Statistical Methods. Vol. 751. John Wiley & Sons, 2013.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#' y <- dietstudy$weightchange[dietstudy$diet == 'Low Fat']
#'
#' lc.signtest.conf <- signtest.conf(x)
#' lf.signtest.conf <- signtest.conf(y)
#'
#' @export signtest.conf
signtest.conf <- function(x, plot = TRUE, conf.level = 0.95){
  # DMD: NOTE: Not currently set up for handling tied values.
  #
  # NOTE: This is taking a success as a negative score,
  # whereas it is more standard to take a success as a
  # positive score.

  n <- length(x)
  x.sort <- sort(x)

  x.median <- if(n %% 2 == 0){
    (x.sort[n/2] + x.sort[n/2 + 1])/2
  }else{
    x.sort[ceiling(n/2)]
  }

  x.min <- x.sort[1]
  x.max <- x.sort[n]

  # See page 66 of Hollander, Wolfe, and Chicken for P-value
  # calculation.

  pconf <- function(mu, lower.tail = TRUE){
    c <- sum(x <= mu)

    if (lower.tail){
      pval <- pbinom(c, n, 0.5) # For right-sided P-values
    }else{
      pval <- pbinom(c-1, n, 0.5, lower.tail = FALSE) # For left-sided P-values
    }

    return(pval)
  }

  pconf <- Vectorize(pconf)

  cconf <- function(mu) {
    c <- sum(x <= mu)

    if (c <= n/2){
      return(1 - 2*pbinom(c-1, n, 1/2))
    }else{
      return(1 - 2*pbinom(n-c, n, 1/2))
    }
  }

  cconf <- Vectorize(cconf)

  # Using procedure from page 80 of Hollander, Wolfe, and Chicken:

  qconf <- function(p) {
    # Using to make qconf give the standard two-sided interval for the sign test.
    if (p >= 0.5){
      x.sort[qbinom(1-p, size = n, prob = 1/2, lower.tail = FALSE) + 1]
    }else{
      x.sort[n + 1 - (qbinom(p, size = n, prob = 1/2, lower.tail = FALSE) + 1)]
    }
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(mu) 1 - cconf(mu)

  scurve <- function(mu) -log2(pcurve(mu))

  out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    display.cconf(out, conf.level = conf.level, xlab = 'Median')

    cat(sprintf("Actual confidence level: %g\n\n", cconf(qconf((1 - conf.level)/2))))
  }

  return(out)
}

#' Confidence Functions for the (Pseudo)Median or Shift via Wilcoxon Tests
#'
#' Confidence functions for the (pseudo)median of a continuous distribution
#' via the Wilcoxon signed-rank test (for one sample) or shift between
#' two continuous distributions via the Wilcoxon rank-sum test (for two samples).
#'
#' @param x a vector of observations for the first sample
#' @param y a vector of observations for the second sample, if applicable
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the median, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Myles Hollander, Douglas A. Wolfe, and Eric Chicken. Nonparametric Statistical Methods. Vol. 751. John Wiley & Sons, 2013.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#' y <- dietstudy$weightchange[dietstudy$diet == 'Low Fat']
#'
#' lc.wc.conf <- wilcox.conf(x)
#' lf.wc.conf <- wilcox.conf(y)
#'
#' comp.wc.conf <- wilcox.conf(x, y)
#'
#' @export wilcox.conf
wilcox.conf <- function(x, y = NULL, plot = TRUE, conf.level = 0.95){
  if (is.null(y)){
    pconf <- Vectorize(function(mu) suppressWarnings(wilcox.test(x, mu = mu, alternative = 'greater')$p.value))

    cconf <- function(mu) abs(2*pconf(mu) - 1)

    # Use wilcox.test()'s upper confidence bound to get quantile function.
    #qconf <- Vectorize(function(p) suppressWarnings(wilcox.test(x, alternative = 'less', conf.level = p, conf.int = TRUE)$conf.int[2]))

    x.range <- range(x)
    x.min <- x.range[1]
    x.max <- x.range[2]

    qconf <- function(p) {
      if (p < pconf(x.min)){
        return(-Inf)
      } else if (p > pconf(x.max)){
        return(Inf)
      } else{
        f.root <- function(x) pconf(x) - p

        q <- uniroot(f.root, interval = c(x.min-1e-10, x.max+1e-10))$root

        return(q)
      }
    }

    qconf <- Vectorize(qconf)

    pcurve <- function(mu) 1 - cconf(mu)

    scurve <- function(mu) -log2(pcurve(mu))

    out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

    if (plot){
      display.cconf(out, conf.level = conf.level, xlab = '(Psuedo)Median')
    }
  }  else{
      pconf <- Vectorize(function(mu) suppressWarnings(wilcox.test(x, y, mu = mu, alternative = 'greater')$p.value))

      cconf <- function(mu) abs(2*pconf(mu) - 1)

      qconf <- Vectorize(function(p) suppressWarnings(wilcox.test(x, y, alternative = 'less', conf.level = p, conf.int = TRUE)$conf.int[2]))

      pcurve <- function(mu) 1 - cconf(mu)

      scurve <- function(mu) -log2(pcurve(mu))

      out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

      if (plot){
        display.cconf(out, conf.level = conf.level, xlab = 'Shift')
      }
  }

  return(out)
}
