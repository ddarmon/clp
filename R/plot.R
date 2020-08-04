#' @export
plot.cconf <- function(obj, conf.level = 0.95, xlim = NULL, xlab = 'Parameter', log = ''){
  if (is.null(xlim)){
    xlim <- obj$qconf(c(0.001, 0.999))
  }

  alpha <- 1-conf.level

  ci <- obj$qconf(c(alpha/2, 1 - alpha/2))

  curve(obj$cconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Curve', n = 2001, lwd = 3, ylim = c(0, 1), log = log)
  segments(x0 = ci[1], x1 = ci[2], y0 = 0.95, col = 'red')
  segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
}

#' @export
plot.dconf <- function(obj, xlim = NULL, xlab = 'Parameter'){
  if (is.null(xlim)){
    xlim <- obj$qconf(c(0.001, 0.999))
  }

  curve(obj$dconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Density', n = 2001, lwd = 3)
}
