#' A Wrapper Function for Plotting Confidence Curves
#'
#' A wrapper for plotting confidence curves using confidence objects
#' returned by *.conf() functions from clp.
#'
#'
#' @param obj a list of confidence functions returned by a *.conf() function
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#' @param xlim the horizontal limits for the plot. If NULL, taken as the larger of
#'             a conf.level confidence interval or a 1 - (1 - conf.level)/2 confidence interval.
#' @param xlab the label for the horizontal axis
#' @param log switch for plotting the horizontal ('x') or vertical ('x') or both ('xy') axes with a
#'            log-scaling
#' @param col the color to plot the confidence curve
#' @param n.points the number of points to plot for the confidence curve
#'
#' @return a plot of the confidence curve from obj
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#'
#' out <- t.conf(x, plot = FALSE)
#'
#' plot.cconf(out, conf.level = 0.999, xlab = 'Average Weight Loss in Low Carb Group (lbs)')
#'
#' @export plot.cconf
plot.cconf <- function(obj, conf.level = 0.95, xlim = NULL, xlab = NULL, log = '', col = 'black', n.points = 2001){
  alpha <- 1-conf.level

  if (is.null(xlim)) {
    if (alpha/2 < 0.001){
      tail.prob <- alpha/4
    }else{
      tail.prob <- 0.001
    }

    xlim <- obj$qconf(c(tail.prob, 1-tail.prob))

    while (xlim[1] == -Inf || xlim[2] == Inf){ # Deal with atoms a +/- Inf
      tail.prob <- tail.prob*1.1

      xlim <- obj$qconf(c(tail.prob, 1-tail.prob))
    }
  }

  if (is.null(xlab)) xlab <- 'Parameter'

  ci <- obj$qconf(c(alpha/2, 1 - alpha/2))

  if (ci[1] == -Inf){
    ci[1] <- xlim[1] - 10
  }
  if (ci[2] == Inf){
    ci[2] <- xlim[2] + 10
  }

  curve(obj$cconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Curve', n = n.points, lwd = 3, ylim = c(0, 1), log = log, col = col)
  segments(x0 = ci[1], x1 = ci[2], y0 = conf.level, col = 'red')
  segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
}

#' A Wrapper Function for Plotting P-Curves
#'
#' A wrapper for plotting P-curves using confidence objects
#' returned by *.conf() functions from clp.
#'
#'
#' @param obj a list of confidence functions returned by a *.conf() function
#' @param conf.level the confidence level for the confidence interval indicated on the P-curve
#' @param xlim the horizontal limits for the plot. If NULL, taken as the larger of
#'             a conf.level confidence interval or a 1 - (1 - conf.level)/2 confidence interval.
#' @param xlab the label for the horizontal axis
#' @param log switch for plotting the horizontal ('x') or vertical ('x') or both ('xy') axes with a
#'            log-scaling
#' @param col the color to plot the confidence curve
#' @param n.points the number of points to plot for the confidence curve
#'
#' @return a plot of the P-curve from obj
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#'
#' out <- t.conf(x, plot = FALSE)
#'
#' plot.pcurve(out, conf.level = 0.999, xlab = 'Average Weight Loss in Low Carb Group (lbs)')
#' @export plot.pcurve
plot.pcurve <- function(obj, conf.level = 0.95, xlim = NULL, xlab = NULL, log = '', col = 'black', n.points = 2001){
  alpha <- 1-conf.level

  if (is.null(xlim)) {
    if (alpha/2 < 0.001){
      tail.prob <- alpha/4
    }else{
      tail.prob <- 0.001
    }

    xlim <- obj$qconf(c(tail.prob, 1-tail.prob))

    while (xlim[1] == -Inf || xlim[2] == Inf){ # Deal with atoms a +/- Inf
      tail.prob <- tail.prob*1.1

      xlim <- obj$qconf(c(tail.prob, 1-tail.prob))
    }
  }

  if (is.null(xlab)) xlab <- 'Parameter'

  ci <- obj$qconf(c(alpha/2, 1 - alpha/2))

  if (ci[1] == -Inf){
    ci[1] <- xlim[1] - 10
  }
  if (ci[2] == Inf){
    ci[2] <- xlim[2] + 10
  }

  curve(obj$pcurve(x), xlim = xlim, xlab = xlab, ylab = expression(italic(P)*'-curve'), n = n.points, lwd = 3, ylim = c(0, 1), log = log, col = col)
  segments(x0 = ci[1], x1 = ci[2], y0 = 1-conf.level, col = 'red')
  segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
}

#' A Wrapper Function for Plotting S-Curves
#'
#' A wrapper for plotting S-curves using confidence objects
#' returned by *.conf() functions from clp.
#'
#'
#' @param obj a list of confidence functions returned by a *.conf() function
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#' @param xlim the horizontal limits for the plot. If NULL, taken as the larger of
#'             a conf.level confidence interval or a 1 - (1 - conf.level)/2 confidence interval.
#' @param xlab the label for the horizontal axis
#' @param log switch for plotting the horizontal ('x') or vertical ('x') or both ('xy') axes with a
#'            log-scaling
#' @param col the color to plot the confidence curve
#' @param n.points the number of points to plot for the confidence curve
#'
#' @return a plot of the S-curve from obj
#'
#' @references  Zad R. Chow and Sander Greenland. "Semantic and Cognitive Tools to Aid Statistical Inference: Replace Confidence and Significance by Compatibility and Surprise." arXiv preprint arXiv:1909.08579 (2019).
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#'
#' out <- t.conf(x, plot = FALSE)
#'
#' plot.scurve(out, conf.level = 0.999, xlab = 'Average Weight Loss in Low Carb Group (lbs)')
#'
#' @export plot.scurve
plot.scurve <- function(obj, conf.level = 0.95, xlim = NULL, xlab = NULL, log = '', col = 'black', n.points = 2001){
  alpha <- 1-conf.level

  if (is.null(xlim)) {
    if (alpha/2 < 0.001){
      tail.prob <- alpha/4
    }else{
      tail.prob <- 0.001
    }

    xlim <- obj$qconf(c(tail.prob, 1-tail.prob))

    while (xlim[1] == -Inf || xlim[2] == Inf){ # Deal with atoms a +/- Inf
      tail.prob <- tail.prob*1.1

      xlim <- obj$qconf(c(tail.prob, 1-tail.prob))
    }
  }

  if (is.null(xlab)) xlab <- 'Parameter'

  ci <- obj$qconf(c(alpha/2, 1 - alpha/2))

  if (ci[1] == -Inf){
    ci[1] <- xlim[1] - 10
  }
  if (ci[2] == Inf){
    ci[2] <- xlim[2] + 10
  }

  curve(obj$scurve(x), xlim = xlim, xlab = xlab, ylab = expression(italic(S)*'-curve (bits)'), n = n.points, lwd = 3, log = log, col = col)
  segments(x0 = ci[1], x1 = ci[2], y0 = -log2(1-conf.level), col = 'red')
  segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
}

#' A Wrapper Function for Plotting Confidence Densities
#'
#' A wrapper for plotting confidence densities using confidence objects
#' returned by *.conf() functions from clp.
#'
#'
#' @param obj a list of confidence functions returned by a *.conf() function
#' @param xlim the horizontal limits for the plot. If NULL, taken as a 99.8\%
#'             confidence interval.
#' @param xlab the label for the horizontal axis
#' @param col the color to plot the confidence curve
#' @param n.points the number of points to plot for the confidence curve
#'
#' @return a plot of the confidence curve from obj
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#'
#' out <- t.conf(x, plot = FALSE)
#'
#' plot.dconf(out, xlab = 'Average Weight Loss in Low Carb Group (lbs)')
#' @export plot.dconf
plot.dconf <- function(obj, xlim = NULL, xlab = NULL, col = 'black', n.points = 2001){
  if (is.null(xlim)) xlim <- obj$qconf(c(0.001, 0.999))
  if (is.null(xlab)) xlab <- 'Parameter'

  curve(obj$dconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Density', n = n.points, lwd = 3, col = col)
}

#' A Wrapper Function for Plotting Confidence Distributions
#'
#' A wrapper for plotting confidence distributions using confidence objects
#' returned by *.conf() functions from clp.
#'
#'
#' @param obj a list of confidence functions returned by a *.conf() function
#' @param compl whether to plot the complementary
#'                      confidence distribution (TRUE) or
#'                      confidence distribution (FALSE, default).
#' @param xlim the horizontal limits for the plot. If NULL, taken as a 99.8\%
#'             confidence interval.
#' @param xlab the label for the horizontal axis
#' @param col the color to plot the confidence curve
#' @param n.points the number of points to plot for the confidence curve
#'
#' @return a plot of the confidence curve from obj
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' data(dietstudy)
#'
#' x <- dietstudy$weightchange[dietstudy$diet == 'Low Carb']
#'
#' out <- t.conf(x, plot = FALSE)
#'
#' plot.pconf(out, xlab = 'Average Weight Loss in Low Carb Group (lbs)')
#'
#' @export plot.pconf
plot.pconf <- function(obj, compl = FALSE, xlim = NULL, xlab = NULL, col = 'black', n.points = 2001){
  if (is.null(xlim)) xlim <- obj$qconf(c(0.001, 0.999))
  if (is.null(xlab)) xlab <- 'Parameter'

  if (compl){
    curve(1 - obj$pconf(x), xlim = xlim, xlab = xlab, ylab = 'Complementary Confidence Distribution', n = n.points, lwd = 3, col = col)
  }else{
    curve(obj$pconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Distribution', n = n.points, lwd = 3, col = col)
  }
}
