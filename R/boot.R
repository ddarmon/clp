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
#' @param formula for use with interfacing to lm, glm, etc.
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
#' # Bootstrap confidence functions for a single mean.
#'
#' t.one.sample <- function(data, id = 1:length(data), ...){
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
#' # Reproduce BCa confidence density from Figure 11.7
#' # of *Computer Age Statistical Inference*.
#'
#' scor <-
#' read.table('https://web.stanford.edu/~hastie/CASI_files/DATA/student_score.txt',
#' header = TRUE)
#'
#' statistic <- function(data, id = 1:nrow(data)){
#'   dat <- data[id, ]
#'
#'   Sigma <- cov(dat)*((nrow(dat)-1)/nrow(dat))
#'
#'   lams <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#'
#'   return(lams[1])
#' }
#'
#' ran.gen <- function(data, mle){
#'   mvrnorm(n = nrow(data), mle$mu, mle$Sigma)
#' }
#'
#' bc <- bcaboot(data = scor, statistic = statistic, B = 8000, sim = "parametric", ran.gen = ran.gen,
#'               mle = list(mu = colMeans(scor),
#'                          Sigma = cov(scor)*((nrow(scor)-1)/nrow(scor))))
#'
#' @export
bcaboot <- function(data, statistic, B = 2000, sim = "ordinary", stratified = FALSE, ran.gen = function(d, p) d, mle = NULL, formula = NULL){
  if (stratified){
    strata <- data[, ncol(data)]

    if (is.null(formula)){
      boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B, sim = sim, ran.gen = ran.gen, mle = mle, parallel = 'multicore', ncpus = parallel::detectCores())
    }else{
      boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B, sim = sim, ran.gen = ran.gen, mle = mle, formula = formula, parallel = 'multicore', ncpus = parallel::detectCores())
    }

  }else{
    if (is.null(formula)){
      boot.out <- boot(data = data, statistic = statistic, R = B, sim = sim, ran.gen = ran.gen, mle = mle, parallel = 'multicore', ncpus = parallel::detectCores())
    }else{
      boot.out <- boot(data = data, statistic = statistic, R = B, sim = sim, ran.gen = ran.gen, mle = mle, formula = formula, parallel = 'multicore', ncpus = parallel::detectCores())
    }
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

  if (is.null(formula)){
    for (i in seq_len(n)) {
      u[i, ] <- statistic(data[-i, ], seq_len(n-1))
    }
  }else{
    for (i in seq_len(n)) {
      u[i, ] <- statistic(data[-i, ], seq_len(n-1), formula = formula)
    }
  }

  t. <- sweep(-u, 2, colMeans(u), "+") * (n - 1)
  a <- (1 / 6) * colSums(t.^3) / (colSums(t.^2))^1.5

  Gn <- list()

  for (i in 1:p){
    Gn[[i]] = stats::ecdf(boot.out$t[, i])
  }

  return(list(t0 = boot.out$t0, t = boot.out$t, Gn = Gn, z0 = z0, a = a))
}

#' Perform Percentile Bootstrap for a Given Statistic
#'
#' Approximate the bootstrap distribution for the percentile bootstrap, given a statistic as would be supplied to
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
#' @param formula for use with interfacing to lm, glm, etc.
#'
#' @return A list containing the sample estimate, bootstrap estimates,
#'         and bootstrap CDF for the percentile bootstrap.
#'
#' @references  Thomas J. DiCiccio and Bradley Efron. "Bootstrap confidence intervals." Statistical Science (1996): 189-212.
#'
#' @examples
#' # Bootstrap confidence functions for a single mean.
#'
#' t.one.sample <- function(data, id = 1:length(data), ...){
#'   dat <- data[id]
#'
#'   d <- mean(dat)
#'
#'   return(d)
#' }
#'
#' data(dietstudy)
#'
#' percboot(data = dietstudy$weightchange[dietstudy$diet == 'Low Carb'],
#'         statistic = t.one.sample,
#'         B = 2000)
#'
#'
#' @export percboot
percboot <- function(data, statistic, B = 2000, sim = "ordinary", stratified = FALSE, ran.gen = function(d, p) d, mle = NULL, formula = NULL){
  if (stratified){
    strata <- data[, ncol(data)]

    if (is.null(formula)){
      boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B, sim = sim, ran.gen = ran.gen, mle = mle, parallel = 'multicore', ncpus = parallel::detectCores())
    }else{
      boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B, sim = sim, ran.gen = ran.gen, mle = mle, formula = formula, parallel = 'multicore', ncpus = parallel::detectCores())
    }

  }else{
    if (is.null(formula)){
      boot.out <- boot(data = data, statistic = statistic, R = B, sim = sim, ran.gen = ran.gen, mle = mle, parallel = 'multicore', ncpus = parallel::detectCores())
    }else{
      boot.out <- boot(data = data, statistic = statistic, R = B, sim = sim, ran.gen = ran.gen, mle = mle, formula = formula, parallel = 'multicore', ncpus = parallel::detectCores())
    }
  }

  p <- length(boot.out$t0) # Dimension of the parameter vector

  Gn <- list()

  for (i in 1:p){
    Gn[[i]] = stats::ecdf(boot.out$t[, i])
  }

  return(list(t0 = boot.out$t0, t = boot.out$t, Gn = Gn))
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

  # Handle when Gn(theta) is 0 or 1.

  which.pinf <- which(Phi.invs == Inf)

  if (length(which.pinf) > 0){
    Hn[which.pinf] <- 1
    warning("Warning: Evaluating BCa confidence distribution at a parameter value greater than the largest bootstrapped estimate.\nInferential statistics may be unreliable.\n")
  }

  which.ninf <- which(Phi.invs == -Inf)

  if (length(which.ninf) > 0){
    Hn[which.ninf] <- 0
    warning("Warning: Evaluating BCa confidence distribution at a parameter value less than the smallest bootstrapped estimate.\nInferential statistics may be unreliable.\n")
  }

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

#' Construct the Confidence Distribution from a Percentile Bootstrap for a Given Statistic
#'
#' Construct the percentile confidence distribution from a bootstrap sample constructed using
#' bcaboot.
#'
#' @param bc an object returned by percboot
#' @param theta the parameter value at which to evaluate the percentile confidence distribution
#' @param param which parameter value to construct the percentile confidence distribution for
#'
#' @return The percentile confidence distribution for param evaluated at theta.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confdist.perc = function(bc, theta, param){
  Gn = bc$Gn[[param]]

  return(Gn(theta))
}

#' Construct the Confidence Density from a Percentile Bootstrap for a Given Statistic
#'
#' Construct the percentile confidence density from a bootstrap sample constructed using
#' percboot.
#'
#' @param bc an object returned by percboot
#' @param param which parameter value to construct the percentile confidence density for
#'
#' @return A spline approximation the percentile confidence density.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confdens.perc = function(bc, param){
  density.out <- density(bc$t[, param], bw = "SJ", n = 1024) # Seems to undersmooth

  gn.perc <- density.out$y
  thetas <- density.out$x

  gn.approx <- approxfun(thetas, gn.perc)

  dconf <- function(thetas){
    g <- gn.approx(thetas)

    g[is.na(g)] <- 0

    return(g)
  }

  return(dconf)
}

#' Construct the Confidence Quantile Function from a Percentile Bootstrap for a Given Statistic
#'
#' Construct the percentile confidence quantile function from a bootstrap sample constructed using
#' bcaboot.
#'
#' @param bc an object returned by percboot
#' @param p the desired probability value at which to evaluate the confidence quantile function
#' @param param which parameter value to construct the percentile confidence quantile function for
#'
#' @return The percentile confidence quantile function for param evaluated at p.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
confquant.perc <- function(bc, p, param){
  num.na <- sum(is.na(bc$t[, param]))

  if (num.na > 0) {
    warning(sprintf('%g of %g bootstrap estimates are NA. Removing them when approximating estimating the bootstrap distribution.', num.na, nrow(bc$t)))
    Qn <- quantile(bc$t[, param], p, na.rm = TRUE)
  }else{
    Qn <- quantile(bc$t[, param], p)
  }


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
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A single value of the difference of the two bootstrapped sample means.
#'
t.two.sample <- function(data, id = 1:nrow(data), ...){
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
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A single value of the bootstrapped sample mean.
#'
t.one.sample <- function(data, id = 1:length(data), ...){
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

  out <- conffuns.from.bcaboot(bc)

  if (plot){
    plot.dconf(out, xlab = xlab)
    plot.cconf(out, conf.level = conf.level, xlab = xlab)
  }

  return(out)
}

lm.beta.for.boot <- function(data, id = 1:nrow(data), formula){
  dat <- data[id, ]

  mod <- lm(formula, data = dat)

  return(mod$coefficients)
}

#' Bootstrapped Confidence Functions for Coefficients of a Linear Model
#'
#' Bootstrapped confidence functions for coefficients of a linear model
#' using BCa and the case resampling bootstrap.
#'
#' @param formula the model to be fitted
#' @param data the data frame containing the data for fitting
#' @param B the number of bootstrap samples used to approximate the
#'          bootstrap distribution
#'
#' @return A list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each coefficient of the linear model, as well as
#'         the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
#' @examples
#' data(fat)
#'
#' beta.boot.conf <- lm.beta.boot.conf(body.fat ~ age + weight + height, data = fat, B = 2000)
#'
#' plot.cconf(beta.boot.conf$weight)
#'
#' @export lm.beta.boot.conf
lm.beta.boot.conf <- function(formula, data, B = 2000){
  bc <- bcaboot(data = data, statistic = lm.beta.for.boot, B = B, formula = formula)

  out <- conffuns.from.bcaboot(bc)

  return(out)
}

lm.sigma.for.boot <- function(data, id = 1:nrow(data), formula){
  dat <- data[id, ]

  mod <- lm(formula, data = dat)

  df <- mod$df.residual

  sig <- sqrt(df*(summary(mod)$sigma)^2/nrow(mod$model))

  return(sig)
}

lm.sigma.boot.conf <- function(formula, data, B = 2000, plot = TRUE, conf.level = 0.95){
  bc <- bcaboot(data = data, statistic = lm.sigma.for.boot, B = B, formula = formula)

  out <- conffuns.from.bcaboot(bc)

  if (plot){
    plot.dconf(out, xlab = 'Noise Standard Deviation')
    plot.cconf(out, conf.level = conf.level, xlab = 'Noise Standard Deviation')
  }

  return(out)
}

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

conffuns.from.percboot.single <- function(bc, ind){
  ind

  pconf <- function(x) confdist.perc(bc, x, ind)
  cconf <- function(x) abs(2*pconf(x) - 1)

  # Only valid within the approximate range of the data!
  dconf <- function(x) confdens.perc(bc, ind)(x)

  qconf <- function(p) confquant.perc(bc, p, ind)

  pcurve <- function(x) 1 - cconf(x)
  scurve <- function(x) -log2(pcurve(x))

  out <- list(pconf = pconf, cconf = cconf, dconf = dconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  return(out)
}

#' Construct Confidence Functions from an Output of bcaboot
#'
#' Construct confidence functions from the output of bcaboot.
#' This is a helper function to make constructing ones own
#' bootstrap confidence functions straightforward.
#'
#' @param bc an output from bcaboot
#'
#' @return A list or list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the parameter associated with the statistic(s) used with
#'         bcaboot.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
#' @examples
#' # Reproduce BCa confidence density from Figure 11.7
#' # of *Computer Age Statistical Inference*.
#'
#' scor <-
#' read.table('https://web.stanford.edu/~hastie/CASI_files/DATA/student_score.txt',
#' header = TRUE)
#'
#' statistic <- function(data, id = 1:nrow(data), ...){
#'   dat <- data[id, ]
#'
#'   Sigma <- cov(dat)*((nrow(dat)-1)/nrow(dat))
#'
#'   lams <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#'
#'   return(lams[1])
#' }
#'
#' ran.gen <- function(data, mle){
#'   mvrnorm(n = nrow(data), mle$mu, mle$Sigma)
#' }
#'
#' bc <- bcaboot(data = scor, statistic = statistic, B = 8000, sim = "parametric", ran.gen = ran.gen,
#'               mle = list(mu = colMeans(scor),
#'                          Sigma = cov(scor)*((nrow(scor)-1)/nrow(scor))))
#'
#' lam.conf <- conffuns.from.bcaboot(bc)
#'
#' plot.cconf(lam.conf, xlab = 'Largest Eigenvalue')
#' plot.dconf(lam.conf, xlab = 'Largest Eigenvalue')
#'
#' lam.conf$qconf(c(0.025, 0.975))
#'
#' @export conffuns.from.bcaboot
conffuns.from.bcaboot <- function(bc){
  if (length(bc$t0) > 1){ # Parameter vector

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

#' Construct Confidence Functions from an Output of percboot
#'
#' Construct confidence functions from the output of percboot.
#' This is a helper function to make constructing ones own
#' bootstrap confidence functions straightforward.
#'
#' @param bc an output from percboot
#'
#' @return A list or list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the parameter associated with the statistic(s) used with percboot.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Bradley Efron and Trevor Hastie. Computer Age Statistical Inference. Vol. 5. Cambridge University Press, 2016.
#'
#' @examples
#' # Reproduce percentile confidence density from Figure 11.7
#' # of *Computer Age Statistical Inference*.
#'
#' scor <-
#' read.table('https://web.stanford.edu/~hastie/CASI_files/DATA/student_score.txt',
#' header = TRUE)
#'
#' statistic <- function(data, id = 1:nrow(data), ...){
#'   dat <- data[id, ]
#'
#'   Sigma <- cov(dat)*((nrow(dat)-1)/nrow(dat))
#'
#'   lams <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#'
#'   return(lams[1])
#' }
#'
#' ran.gen <- function(data, mle){
#'   mvrnorm(n = nrow(data), mle$mu, mle$Sigma)
#' }
#'
#' bc <- percboot(data = scor, statistic = statistic, B = 8000, sim = "parametric", ran.gen = ran.gen,
#'               mle = list(mu = colMeans(scor),
#'                          Sigma = cov(scor)*((nrow(scor)-1)/nrow(scor))))
#'
#' lam.conf <- conffuns.from.percboot(bc)
#'
#' plot.cconf(lam.conf, xlab = 'Largest Eigenvalue')
#' plot.dconf(lam.conf, xlab = 'Largest Eigenvalue')
#'
#' lam.conf$qconf(c(0.025, 0.975))
#'
#' @export conffuns.from.percboot
conffuns.from.percboot <- function(bc){
  if (is.null(bc$Gn)){ # A boot object made directly by the boot package
    Gn <- list()

    for (i in 1:length(bc$t0)){
      Gn[[i]] = stats::ecdf(bc$t[, i])
    }

    bc$Gn <- Gn
  }

  if (length(bc$t0) > 1){ # Parameter vector

    out <- list()

    if (is.null(names(bc$t0))){
      theta.names <- paste0('theta', 1:length(bc$t0))
    }else{
      theta.names <- names(bc$t0)
    }

    for (ind in 1:length(theta.names)){
      out[[theta.names[ind]]] <- conffuns.from.percboot.single(bc, ind)
    }

    return(out)
  }else{ # Single parameter

    out <- conffuns.from.percboot.single(bc, 1)

    return(out)
  }
}

# Function to generate confidence functions via bootstrapping
# adjusted for multiple comparisons using the method of Rudolf Beran (1988):
#
# Rudolf Beran. "Balanced simultaneous confidence sets." Journal of the American Statistical Association 83.403 (1988): 679-686.
make.beran.multicomp.obj <- function(theta.boot, K, bca.params = NULL){
  Gn <- ecdf(theta.boot)

  if (is.null(bca.params)){
    cconf <- function(theta) K(abs(2*Gn(theta) - 1))

    # Median of the confidence distribution is the median of the
    # bootstrap estimates.

    theta.med <- median(theta.boot)
  }else{
    cconf <- function(theta) {
      Phi.invs = qnorm(Gn(theta))

      # The BCa confidence distribution
      Hn = pnorm((Phi.invs - bca.params$z0)/(1 + bca.params$a*(Phi.invs - bca.params$z0)) - bca.params$z0)

      # Handle when Gn(theta) is 0 or 1.

      which.pinf <- which(Phi.invs == Inf)

      if (length(which.pinf) > 0){
        Hn[which.pinf] <- 1
        #warning("Warning: Evaluating BCa confidence distribution at a parameter value greater than the largest bootstrapped estimate.\nInferential statistics may be unreliable.\n")
      }

      which.ninf <- which(Phi.invs == -Inf)

      if (length(which.ninf) > 0){
        Hn[which.ninf] <- 0
        #warning("Warning: Evaluating BCa confidence distribution at a parameter value less than the smallest bootstrapped estimate.\nInferential statistics may be unreliable.\n")
      }
      return(K(abs(2*Hn - 1)))
    }

    # Median of the BCa confidence distribution is
    #   \psi_{1/2} =  \hat{G}^{-1}\left( \Phi\left( \frac{z_{0}(2 - az_{0})}{1 - a z_{0}} \right) \right)

    inner <- pnorm(bca.params$z0*(2-bca.params$a*bca.params$z0)/(1 - bca.params$a*bca.params$z0))

    theta.med <- as.numeric(quantile(theta.boot, probs = inner))
  }

  # Might want to adjust this too, for when z0 is large.

  theta.min <- min(theta.boot) - 1 # Need - 1 since Fhat(theta.min) = 1/n, and not 0.
  theta.max <- max(theta.boot)

  pconf <- function(theta){
    if (theta <= theta.med){
      return(0.5 - 0.5*cconf(theta))
    }else{
      return(0.5 + 0.5*cconf(theta))
    }
  }

  pconf <- Vectorize(pconf)

  qconf <- function(p){
    fun.root <- function(x) pconf(x) - p

    uniroot(fun.root, interval = c(theta.min, theta.max))$root
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(theta) 1 - cconf(theta)

  scurve <- function(theta) -log2(pcurve(theta))

  out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  return(out)
}

#' Simultaneous Confidence Functions via Bootstrapping.
#'
#' Confidence functions via bootstrapping based on the
#' percentile confidence interval, with balancing across
#' parameters with Beran's method for simultaneous confidence sets.
#'
#' @param bc an object returned by either bcaboot or percboot
#' @param which the indices for the parameters to include
#' @param interval.type which interval to use in constructing the confidence curve, one of 'percentile' or 'bca'
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each parameter with balanced correction for multiple comparisons.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, Likelihood, Probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Tore Schweder. "Confidence nets for curves." Advances In Statistical Modeling And Inference: Essays in Honor of Kjell A Doksum. 2007. 593-609.
#'
#'              Rudolf Beran. "Balanced simultaneous confidence sets." Journal of the American Statistical Association 83.403 (1988): 679-686.
#'
#' @examples
#' data(fat)
#'
#' formula <- body.fat ~ age + weight + height
#'
#' statistic <- function(data, id = 1:nrow(data)){
#'  mod <- lm(formula, data = data[id, ])
#'
#'  return(coef(mod))
#' }
#'
#' bc <- bcaboot(fat, statistic, B = 2000)
#'
#' beta.conf.simul <- bootstrap.beran.conf(bc, which = 2:4)
#'
#' for (nam in names(beta.conf.simul)){
#'  plot.cconf(beta.conf.simul[[nam]], xlab = nam)
#' }
#'
#' @export bootstrap.beran.conf
bootstrap.beran.conf <- function(bc, which = NULL, interval.type = 'bca'){
  if (is.null(which)){
    which <- 1:ncol(bc$t)
  }

  theta.hat  <- bc$t0[which]
  theta.boot <- bc$t[, which]

  if (interval.type == 'bca'){
    if (is.null(bc$z0) || is.null(bc$a)){
      warning('The boot object does not have the necessary BCa parameters. Using the percentile interval confidence distribution instead.')
      interval.type <- 'percentile'
    }else{
      z0.sub <- bc$z0[which]
      a.sub <- bc$a[which]
    }
  }

  nam <- names(theta.hat)

  if (is.null(nam)){
    nam <- as.character(which)
  }

  V <- matrix(nrow = nrow(theta.boot), ncol = ncol(theta.boot))

  for (j in 1:ncol(theta.boot)){
    Hn <- ecdf(theta.boot[, j])

    V[, j] <- abs(2*Hn(theta.boot[, j]) - 1)
  }

  Vmax <- apply(V, 1, max)
  K <- ecdf(Vmax)

  out <- list()

  if (interval.type == 'bca'){
    for (j in 1:ncol(theta.boot)){
      out[[nam[j]]] <- make.beran.multicomp.obj(theta.boot[, j], K, list(z0 = z0.sub[j], a = a.sub[j]))
    }
  }else{
    for (j in 1:ncol(theta.boot)){
      out[[nam[j]]] <- make.beran.multicomp.obj(theta.boot[, j], K)
    }
  }

  return(out)
}
