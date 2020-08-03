#' Confidence functions for t-Test
#'
#' Computes confidence functions for one and two sample t-tests on vectors of data.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param plot a logical indicating whether you want the plots to be generated.
#' @return A list containing the confidence functions.
#' @examples
#' t.confdist(x = 1:10, y = 7:20)


t.conf <- function(x, y = NULL, paired = FALSE, plot = TRUE, conf.level = 0.95) {

  if(!is.null(y)){
    if(paired)
      xok <- yok <- complete.cases(x,y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  } else {
    if (paired) stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if (paired) {
    x <- x-y
    y <- NULL
  }

  xbar <- mean(x)
  sx <- sd(x)
  nx <- length(x)
  vx <- var(x)

  if(is.null(y)){
    se.xbar <- sx/sqrt(nx)

    # Confidence Distribution
    pconf <- function(mu) pt((mu - xbar)/se.xbar, df = (nx-1))

    # Confidence Density
    dconf <- function(mu) dt((mu - xbar)/se.xbar, df = nx-1)

    # Confidence Curve
    cconf <- function(mu) abs(2*pconf(mu) - 1)

    # Confidence Quantile
    qconf <- function(p) xbar + se.xbar*qt(p, df = nx - 1)

    # P-curve
    pcurve <- function(mu) 1 - cconf(mu)

    out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve)

    if (plot){
      plot.dconf(out, xlab = 'Mean')
      plot.cconf(out, conf.level = conf.level, xlab = 'Mean')
    }

    return(out)

  } else if (!is.null(y)){

    ybar <- mean(y)
    sy <- sd(y)
    ny <- length(y)
    vy <- var(y)
    stderrx <- sqrt(vx/nx)
    stderry <- sqrt(vy/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df.welch <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))

    diff <- xbar - ybar
    se.diff <- sqrt(sx^2/nx + sy^2/ny)

    # Confidence Distribution
    pconf <- function(delta) pt((delta-diff)/se.diff, df = df.welch)

    # Confidence Density
    dconf <- function(delta) dt((delta-diff)/se.diff, df = df.welch)/se.diff

    # Confidence Curve
    cconf <- function(delta) abs(2*pconf(delta) - 1)

    # Confidence Quantile
    qconf <- function(p) diff + se.diff*qt(p, df = df.welch)

    # P-curve
    pcurve <- function(delta) 1-cconf(delta)

    out <- list(pconf = pconf, dconf = dconf, cconf =  cconf, qconf = qconf, pcurve = pcurve)

    if (plot){
      plot.dconf(out, xlab = 'Mean[1] - Mean[2]')
      plot.cconf(out, conf.level = conf.level, xlab = 'Mean[1] - Mean[2]')
    }

    return(out)

  }
}

#' Confidence functions for t-Test with summary data
#'
#' Computes confidence functions for one and two sample t-tests on summary data.
#'
#' @param mean a (non-empty) numeric vector of data values.
#' @param sd a (non-empty) numeric vector of data values.
#' @param n a (non-empty) numeric vector of data values.
#' @param plot a logical indicating whether you want the plots to be generated.
#' @return A list containing the confidence functions.
#' @examples
#' t.confdist(c(3,5), c(1,.5), c(10, 12))

t.conf.summary <- function(mean, sd, n, plot = TRUE, conf.level = 0.95){

  stopifnot(length(mean)==length(sd) && length(mean)==length(n))

  if (length(mean)==1){
    xbar <- mean[1]
    sx <- sd[1]
    nx <- n[1]
    vx <- sx^2

    se.xbar <- sx/sqrt(nx)

    # Confidence Distribution
    pconf <- function(mu) pt((mu - xbar)/se.xbar, df = (nx-1))

    # Confidence Density
    dconf <- function(mu) dt((mu - xbar)/se.xbar, df = nx-1)

    # Confidence Curve
    cconf <- function(mu) abs(2*pconf(mu) - 1)

    # Confidence Quantile
    qconf <- function(p) xbar + se.xbar*qt(p, df = nx - 1)

    # P-curve
    pcurve <- function(mu) 1 - cconf(mu)

    out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve)

    if (plot){
      plot.dconf(out, xlab = 'Mean')
      plot.cconf(out, conf.level = conf.level, xlab = 'Mean')
    }

    return(out)

  } else if (length(mean)==2){
    xbar <- mean[1]
    sx <- sd[1]
    nx <- n[1]
    vx <- sx^2
    ybar <- mean[2]
    sy <- sd[2]
    ny <- n[2]
    vy <- sy^2
    stderrx <- sqrt(vx/nx)
    stderry <- sqrt(vy/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df.welch <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))

    diff <- xbar - ybar
    se.diff <- sqrt(sx^2/nx + sy^2/ny)

    # Confidence Distribution
    pconf <- function(delta) pt((delta-diff)/se.diff, df = df.welch)

    # Confidence Density
    dconf <- function(delta) dt((delta-diff)/se.diff, df = df.welch)/se.diff

    # Confidence Curve
    cconf <- function(delta) abs(2*pconf(delta) - 1)

    # Confidence Quantile
    qconf <- function(p) diff + se.diff*qt(p, df = df.welch)

    # P-curve
    pcurve <- function(delta) 1-cconf(delta)

    out <- list(pconf = pconf, dconf = dconf, cconf =  cconf, qconf = qconf, pcurve = pcurve)

    if (plot){
      plot.dconf(out, xlab = 'Mean[1] - Mean[2]')
      plot.cconf(out, conf.level = conf.level, xlab = 'Mean[1] - Mean[2]')
    }

    return(out)

  }else{
    stop("Vector lengths must be 1 or 2.")
  }

}
