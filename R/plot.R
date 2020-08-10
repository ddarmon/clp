#' @export plot.cconf
plot.cconf <- function(obj, conf.level = 0.95, xlim = NULL, xlab = NULL, log = '', which = 1, col = 'black'){
  alpha <- 1-conf.level

  if (class(obj) == 'lm.beta.conf'){
    if (is.null(xlim)) {
      if (alpha/2 < 0.001){
        xlim <- obj$qconf(which)(c(alpha/4, 1-alpha/4))
      }else{
        xlim <- obj$qconf(which)(c(0.001, 0.999))
      }
    }
    if (is.null(xlab)) {
      if (is.numeric(which)){
        xlab <- obj$variable[which]
      }else{
        xlab <- which
      }
    }

    ci <- obj$qconf(which)(c(alpha/2, 1 - alpha/2))
    curve((obj$cconf(which))(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Curve', n = 2001, lwd = 3, ylim = c(0, 1), log = log, col = col)
    segments(x0 = ci[1], x1 = ci[2], y0 = conf.level, col = 'red')
    segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
  }else{
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

    curve(obj$cconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Curve', n = 2001, lwd = 3, ylim = c(0, 1), log = log, col = col)
    segments(x0 = ci[1], x1 = ci[2], y0 = conf.level, col = 'red')
    segments(x0 = ci[1], x1 = ci[2], y0 = -0.025, col = 'red', lwd = 4)
  }
}

#' @export plot.dconf
plot.dconf <- function(obj, xlim = NULL, xlab = NULL, which = 1, col = 'black'){
  if (class(obj) == 'lm.beta.conf'){
    if (is.null(xlim)) xlim <- (obj$qconf(which))(c(0.001, 0.999))
    if (is.null(xlab)) {
      if (is.numeric(which)){
        xlab <- obj$variable.names[which]
      }else{
        xlab <- which
      }
    }

    curve(obj$dconf(which)(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Density', n = 2001, lwd = 3, col = col)
  }else{
    if (is.null(xlim)) xlim <- obj$qconf(c(0.001, 0.999))
    if (is.null(xlab)) xlab <- 'Parameter'

    curve(obj$dconf(x), xlim = xlim, xlab = xlab, ylab = 'Confidence Density', n = 2001, lwd = 3, col = col)
  }
}
