#' Perform BCa Bootstrap for a Given Statistic
#'
#' Approximate the bootstrap distribution and estimate bias and acceleration
#' parameters for the BCa bootstrap, given a statistic as would be supplied to
#' boot from the boot package.
#'
#' @param data the dataframe or data matrix
#' @param statistic a function that computes the statistic from data
#' @param B the number of bootstrap replicates
#' @param sim either "ordinary" (for case resampling bootstrap) or "parametric" (for parametric bootstrap)
#' @param stratified whether or not (default) to use the use stratified
#'                   sampling for the bootstrapping.
#' @param ran.gen a function returning a random sample, for the parametric bootstrap
#' @param mle the maximum likelihood estimate from the original sample, for the parametric bootstrap
#'
#' @return A list containing the sample estimate, bootstrap estimates,
#'         bootstrap CDF, bias, and acceleration parameters for the
#'         BCa bootstrap.
#'
#' @references  Bradley Efron. "Better bootstrap confidence intervals." Journal of the American Statistical Association 82.397 (1987): 171-185.
#'
#'              Thomas J. DiCiccio and Bradley Efron. "Bootstrap confidence intervals." Statistical Science (1996): 189-212.
#'
#' @examples
#' t.one.sample <- function(data, id = 1:length(data)){
#'   dat <- data[id]
#'
#'   d <- mean(dat)
#'
#'   return(d)
#' }
#'
#' data(dietstudy)
#'
#' bcaboot(data = dietstudy$weightchange[dietstudy$diet == 'Low Carb'],
#'         statistic = t.one.sample,
#'         B = 2000)
#'
#' @export
bcaboot <- function(data, statistic, B = 2000, sim = "ordinary", stratified = FALSE, ran.gen = function(d, p) d, mle = NULL){
  if (stratified){
    strata <- data[, ncol(data)]

    boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B, sim = sim, ran.gen = ran.gen, mle = mle)
  }else{
    boot.out <- boot(data = data, statistic = statistic, R = B, sim = sim, ran.gen = ran.gen, mle = mle)
  }


  p <- length(boot.out$t0) # Dimension of the parameter vector

  # Compute bias adjustment:
  z0 <- rep(0, p)

  for (i in 1:p){
    # Handle possible NAs:

    comp <- boot.out$t[, i] <= boot.out$t0[i]
    B.wo.na <- sum(!is.na(comp))

    z0[i] <- qnorm(sum(comp, na.rm = TRUE)/(B.wo.na))
  }

  n <- nrow(data)

  if (is.null(n)){
    n <- length(data)
    # Reshape data into a matrix

    data <- matrix(data, nrow = n)
  }

  # Compute acceleration adjustment:
  u <- matrix(rep(0, n*p), nrow = n)

  n1 <- sqrt(n * (n - 1))

  for (i in seq_len(n)) {
    u[i, ] <- statistic(data[-i, ], seq_len(n-1))
  }
  t. <- sweep(-u, 2, colMeans(u), "+") * (n - 1)
  a <- (1 / 6) * colSums(t.^3) / (colSums(t.^2))^1.5

  Gn <- list()

  for (i in 1:p){
    Gn[[i]] = stats::ecdf(boot.out$t[, i])
  }

  return(list(t0 = boot.out$t0, t = boot.out$t, Gn = Gn, z0 = z0, a = a))
}

#' Construct the Confidence Distribution from a BCa Bootstrap for a Given Statistic
#'
#' Construct the BCa confidence distribution from a bootstrap sample constructed using
#' bcaboot.
#'
#' @param bc an object returned by bcaboot
#' @param theta the parameter value at which to evaluate the BCa confidence distribution
#' @param param which parameter value to construct the BCa confidence distribution for
#'
#' @return The BCa confidence distribution for param evaluated at theta.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confdist = function(bc, theta, param){
  Gn = bc$Gn[[param]]
  Phi.invs = qnorm(Gn(theta))

  # The BCa confidence distribution
  Hn = pnorm((Phi.invs - bc$z0[param])/(1 + bc$a[param]*(Phi.invs - bc$z0[param])) - bc$z0[param])

  return(Hn)
}

#' Construct the Confidence Density from a BCa Bootstrap for a Given Statistic
#'
#' Construct the BCa confidence density from a bootstrap sample constructed using
#' bcaboot.
#'
#' @param bc an object returned by bcaboot
#' @param param which parameter value to construct the BCa confidence density for
#'
#' @return A spline approximation the BCa confidence density.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confdens = function(bc, param){
  # See page 202 of *Computer Age Statistical Inference* by
  # Efron and Hastie for the appropriate weights for the
  # BCa confidence density.

  z0 <- bc$z0[param]; a <- bc$a[param]

  Gn <- bc$Gn[[param]]

  w = function(theta, Gn, z0, a){
    ztheta = qnorm(Gn(theta)) - z0

    bca.fac = dnorm(ztheta/(1+a*ztheta) - z0)/((1 + a*ztheta)^2*dnorm(ztheta + z0))

    return(bca.fac)
  }

  Ws <- w(bc$t[, param], Gn, z0, a)

  Ws[is.na(Ws)] <- 0

  Ws <- Ws/sum(Ws)

  density.out <- density(bc$t[, param], bw = "SJ", weights = Ws, n = 1024) # Seems to undersmooth
  # density.out = density(bc$t[, param], bw = "bcv", weights = Ws, n = 1024) # Seems to oversmooth

  gn.bca <- density.out$y
  thetas <- density.out$x

  gn.approx <- approxfun(thetas, gn.bca)

  dconf <- function(thetas){
    g <- gn.approx(thetas)

    g[is.na(g)] <- 0

    return(g)
  }

  return(dconf)
}

#' Construct the Confidence Quantile Function from a BCa Bootstrap for a Given Statistic
#'
#' Construct the BCa confidence quantile function from a bootstrap sample constructed using
#' bcaboot.
#'
#' @param bc an object returned by bcaboot
#' @param p the desired probability value at which to evaluate the confidence quantile function
#' @param param which parameter value to construct the BCa confidence quantile function for
#'
#' @return The BCa confidence quantile function for param evaluated at p.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confquant <- function(bc, p, param){
  Gn <- bc$Gn[[param]]
  z0 <- bc$z0[param]

  a <- bc$a[param]

  z <- qnorm(p)

  zp <- z0 + (z0 + z)/(1 - a*(z0 + z))

  Qn <- quantile(bc$t[, param], pnorm(zp))

  return(Qn)
}

#' Statistic for Bootstrapped Two-sample t-test.
#'
#' The function to pass to bcaboot to perform the two-sample t-test
#' via bootstrapping.
#'
#' @param data (m + n) x 2 matrix containing the first and second samples,
#'             with the first column containing the data values and the second
#'             column containing an indicator for which sample the data values
#'             are from.
#' @param id the id variable used by boot
#'
#' @return A single value of the difference of the two bootstrapped sample means.
#'
t.two.sample <- function(data, id = 1:nrow(data)){
  dat <- data[id, ]

  d <- mean(dat[dat[, 2] == 1, 1]) - mean(dat[dat[, 2] == 2, 1])

  return(d)
}

#' Statistic for Bootstrapped One-sample t-test.
#'
#' The function to pass to bcaboot to perform the one-sample t-test
#' via bootstrapping.
#'
#' @param data a vector containing the sample
#' @param id the id variable used by boot
#'
#' @return A single value of the bootstrapped sample mean.
#'
t.one.sample <- function(data, id = 1:length(data)){
  dat <- data[id]

  d <- mean(dat)

  return(d)
}

#' Bootstrapped Confidence Functions for One or Two Means
#'
#' Bootstrapped confidence functions for a single mean or the difference
#' of two means using the BCa bootstrap.
#'
#' @param x a vector containing the first sample
#' @param y a vector containing the second sample (optional)
#' @param B the number of bootstrap samples used to approximate the
#'          bootstrap distribution
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for a single mean or the difference of two means, as well as
#'         the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
#' @examples
#'
#' data(dietstudy)
#'
#' # One mean
#'
#' t.boot.conf(x = dietstudy$weightchange[dietstudy$diet == 'Low Carb'],
#'             B = 2000)
#'
#' # Two means
#'
#' t.boot.conf(x = dietstudy$weightchange[dietstudy$diet == 'Low Carb'],
#'             y = dietstudy$weightchange[dietstudy$diet == 'Low Fat'],
#'             B = 2000)
#'
#' @export t.boot.conf
t.boot.conf <- function(x, y = NULL, B = 2000, plot = TRUE, conf.level = 0.95){

  if (is.null(y)){
    dat <- x
    bc <- bcaboot(dat, t.one.sample, B = B)

    xlab <- 'mean'
  }else{
    dat <- cbind(c(x, y), c(rep(1, length(x)), rep(2, length(y))))
    bc <- bcaboot(dat, t.two.sample, B = B, stratified = TRUE)

    xlab <- 'mean[1] - mean[2]'
  }

  out <- conffuns.from.bcaboot(bc, 1)

  if (plot){
    plot.dconf(out, xlab = xlab)
    plot.cconf(out, conf.level = conf.level, xlab = xlab)
  }

  return(out)
}

#' @export conffuns.from.bcaboot.single
conffuns.from.bcaboot.single <- function(bc, ind){
  ind

  pconf <- function(x) confdist(bc, x, ind)
  cconf <- function(x) abs(2*pconf(x) - 1)

  # Only valid within the approximate range of the data!
  dconf <- function(x) confdens(bc, ind)(x)

  qconf <- function(p) confquant(bc, p, ind)

  pcurve <- function(x) 1 - cconf(x)
  scurve <- function(x) -log2(pcurve(x))

  out <- list(pconf = pconf, cconf = cconf, dconf = dconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  return(out)
}

#' @export conffuns.from.bcaboot
conffuns.from.bcaboot <- function(bc){
  if (length(bc$t0) > 1){ # Parameter vector
    # DMD: This may not work without first calling items first.

    out <- list()

    if (is.null(names(bc$t0))){
      theta.names <- paste0('theta', 1:length(bc$t0))
    }else{
      theta.names <- names(bc$t0)
    }

    for (ind in 1:length(theta.names)){
      out[[theta.names[ind]]] <- conffuns.from.bcaboot.single(bc, ind)
    }

    return(out)
  }else{ # Single parameter

    out <- conffuns.from.bcaboot.single(bc, 1)

    return(out)
  }
}
