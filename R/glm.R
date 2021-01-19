#' Confidence Functions for the Expected Response of a GLM without Dispersion
#'
#' Confidence functions for the expected response of a GLM without dispersion.
#'
#' @param mod an output from glm
#' @param x a vector containing the predictor values for which the expected
#'          response is desired
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the expected response at x, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @examples
#' # Low birth weight example from page 38 of *Confidence, Likelihood, Probability*,
#'
#' data(lbw)
#'
#' mod <- glm(low ~ weight + age + black + other + smoker, data = lbw, family = binomial)
#'
#' x.black <- c(1, 50, 23.238, 1, 0, 1)
#'
#' pred.black <- glm.lincom.conf(mod, x.black)
#'
#' @importFrom utils tail
#' @export
glm.lincom.conf <- function(mod, x, plot = TRUE, conf.level = 0.95){
  # PDF / PMF from GLM
  dmod <- enrichwith::get_dmodel_function(mod)

  # Derivative of log-likelihood
  score <- enrichwith::get_score_function(mod)

  # Log-likelihood evaluated at beta
  ell.beta <- function(beta){
    as.numeric(sum(log(dmod(coefficients = beta))))
  }

  beta.mle <- coef(mod)
  ell.mle <- ell.beta(beta.mle)

  m.mle <- sum(x*beta.mle)

  # Get out feasible direction for the
  # constraint x*beta = m via QR decomp.

  qr.out <- pracma::householder(t(t(x)))

  Fcon <- qr.out$Q[, 2:length(x)]

  # Log-likelihood in terms of the feasible
  # direction.

  ell.con <- function(beta, xhat){
    beta <- matrix(beta, nrow = length(beta))

    arg <- Fcon%*%beta + xhat
    return(ell.beta(arg))
  }

  # Score in terms of the feasible direction.

  score.con <- function(beta, xhat){
    beta <- matrix(beta, nrow = length(beta))

    arg <- Fcon%*%beta + xhat

    t(Fcon)%*%score(coefficients = arg)
  }

  # The profile likelihood function.

  profile.lik <- function(m, beta0){
    opt.out <- optim(beta0, ell.con, score.con, method = 'BFGS', xhat = xhat, control = list(fnscale = -1))
    value <- opt.out$value
    beta.con <- opt.out$par

    return(list(value = value, beta.con = beta.con))
  }

  # Use procedure as described on page 23 of
  #
  # *Fitting Linear Mixed-Effects Models Using lme4*
  # (Bates, Machler, Bolker, and Walker)
  #
  # to perform the profiling.

  ms <- c(m.mle)
  zs <- c(0)

  z.max <- qnorm(0.9999)
  D.z <- z.max / 16

  # Probably should take advantage of asymptotic results here
  # rather than use a fudge factor.

  new.m <- if (ms == 0) 0.001 else 1.01*ms

  ms <- c(ms, new.m)

  # Get value at initial new point:

  xhat <- qr.out$Q[, 1]*(ms[2]/qr.out$R[1])

  beta0 <- t(Fcon)%*%(beta.mle - xhat)

  prof.out <- profile.lik(ms[2], beta0)

  pll <- prof.out$value
  beta0 <- prof.out$beta.con

  zs <- c(zs, sqrt(2*(ell.mle - pll)))

  # Profile in positive direction.

  i <- 2
  while (tail(zs, 1) < z.max){
    DzDpsi <- (zs[i] - zs[i - 1])/(ms[i] - ms[i-1])

    dpsi <- D.z / DzDpsi

    ms <- c(ms, tail(ms, 1) + dpsi)

    xhat <- qr.out$Q[, 1]*(ms[length(ms)]/qr.out$R[1])

    prof.out <- profile.lik(ms[length(ms)], beta0)

    pll <- prof.out$value
    theta0 <- prof.out$beta.con

    zs <- c(zs, sqrt(2*(ell.mle - pll)))

    i <- i + 1
  }

  # Profile in negative direction.

  i <- 2
  while (zs[1] > -z.max){
    DzDpsi <- (zs[i-1] - zs[i])/(ms[i-1] - ms[i])

    dpsi <- D.z / DzDpsi

    ms <- c(ms[1] - dpsi, ms)

    xhat <- qr.out$Q[, 1]*(ms[1]/qr.out$R[1])

    prof.out <- profile.lik(ms[1], beta0)

    pll <- prof.out$value
    theta0 <- prof.out$beta.con

    zs <- c(-sqrt(2*(ell.mle - pll)), zs)

    i <- i + 1
  }

  pll <- ell.mle - (zs^2)/2

  profile.lik.fun <- splinefun(mod$family$linkinv(ms), pll)

  dev.fun <- function(m){
    2*(ell.beta(beta.mle) - profile.lik.fun(m))
  }

  resp.vals <- mod$family$linkinv(ms)

  # For direct comparison with figures from CLP:
  # par(mfrow = c(1, 2))
  # curve(profile.lik.fun, from = min(resp.vals), to = max(resp.vals),
  #       xlab = 'Expected Response', ylab = 'Profile Likelihood')
  # curve(dev.fun,  from = min(resp.vals), to = max(resp.vals),
  #       xlab = 'Expected Response', ylab = 'Deviance')
  # abline(h = qchisq(c(0.5, 0.9, 0.95, 0.99), 1), lty = 3)

  signed.sqrt.dev.fun <- function(m) sign(m - mod$family$linkinv(m.mle))*sqrt(dev.fun(m))

  cconf <- function(m) pchisq(dev.fun(m), 1)
  pconf <- function(m) pnorm(signed.sqrt.dev.fun(m))

  dconf <-  function(m, dx = 1e-5) (pconf(m + dx) - pconf(m - dx))/(2*dx)

  qconf <- function(p) {
    fun.root <- function(z) pconf(z) - p

    rang.ms <- range(mod$family$linkinv(ms))

    return(uniroot(fun.root, interval = c(rang.ms[1], rang.ms[2]))$root) # Might need better upperbound here!
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(m) 1 - cconf(m)

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, profile.lik = profile.lik.fun, deviance = dev.fun)

  if (plot){
    plot.dconf(out, xlab = 'Expected Response')
    plot.cconf(out, conf.level = conf.level, xlab = 'Expected Response')
  }

  return(out)
}

# DMD: Probably should consolidate this into one "master function"
# that chooses the method to use based on whether or not the
# model has a dispersion parameter.

#' Confidence Functions for the Expected Response of a GLM with Dispersion
#'
#' Confidence functions for the expected response of a GLM with dispersion.
#'
#' @param mod an output from glm
#' @param x a vector containing the predictor values for which the expected
#'          response is desired
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the expected response at x, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @export
glm.lincom.conf.disp <- function(mod, x, plot = TRUE, conf.level = 0.95){
  n <- nrow(mod$model)
  p <- length(mod$coefficients)

  # PDF / PMF from GLM
  dmod <- enrichwith::get_dmodel_function(mod)

  # Derivative of log-likelihood
  score <- enrichwith::get_score_function(mod)

  # Log-likelihood evaluated at theta
  ell.theta <- function(theta){
    as.numeric(sum(log(dmod(coefficients = theta[1:(length(theta)-1)], dispersion = theta[length(theta)]))))
  }

  score.theta <- function(theta){
    as.numeric(score(coefficients = theta[1:(length(theta)-1)], dispersion = theta[length(theta)]))
  }

  # MLE of Dispersion Parameter
  phi.hat <- enrichwith::enrich(mod, with = "mle of dispersion")$dispersion_mle

  theta.mle <- c(coef(mod), phi.hat)

  ell.mle <- ell.theta(theta.mle)

  m.mle <- sum(x*theta.mle[1:length(x)])

  dev.mle <- deviance(mod)

  # Get out feasible direction for the
  # constraint x*theta = m via QR decomp.

  x.aug <- c(x, 0)

  qr.out <- pracma::householder(t(t(x.aug)))

  Fcon <- qr.out$Q[, 2:length(x.aug)]

  # Log-likelihood in terms of the feasible
  # direction.

  ell.con <- function(theta, xhat){
    theta <- matrix(theta, nrow = length(theta))

    arg <- Fcon%*%theta + xhat
    arg[length(arg)] <- arg[length(arg)]^2 # Square to enforce positivity

    return(ell.theta(arg))
  }

  # Score in terms of the feasible direction.

  score.con <- function(theta, xhat){
    theta <- matrix(theta, nrow = length(theta))

    arg <- Fcon%*%theta + xhat
    arg[length(arg)] <- arg[length(arg)]^2 # Square to enforce positivity

    t(Fcon)%*%score(coefficients = arg[1:(length(arg)-1)], dispersion = arg[length(arg)])
  }

  # The profile likelihood function.

  profile.lik <- function(m, theta0){
    # Suppressing warning for NAs that can occur for gamma GLMs.
    suppressWarnings(opt.out <- optim(theta0, ell.con, score.con, method = 'BFGS', xhat = xhat, control = list(fnscale = -1)))
    value <- opt.out$value
    theta.con <- opt.out$par

    return(list(value = value, theta.con = theta.con))
  }

  # Use procedure as described on page 23 of
  #
  # *Fitting Linear Mixed-Effects Models Using lme4*
  # (Bates, Machler, Bolker, and Walker)
  #
  # to perform the profiling.

  ms <- c(m.mle)
  zs <- c(0)

  z.max <- qt(0.9999, n - p)
  D.z <- z.max / 16

  # Probably should take advantage of asymptotic results here
  # rather than use a fudge factor.

  new.m <- if (ms == 0) 0.001 else 1.01*ms

  ms <- c(ms, new.m)

  # Get values at initial new point:

  xhat <- qr.out$Q[, 1]*(ms[2]/qr.out$R[1])

  # Do these two lines in non-disp. code as well.

  theta0 <- t(Fcon)%*%(theta.mle - xhat)

  theta0[length(theta0)] <- sqrt(theta0[length(theta0)])

  # theta0 <- theta.mle[2:length(theta.mle)]

  prof.out <- profile.lik(ms[2], theta0)

  pll <- prof.out$value
  theta0 <- prof.out$theta.con

  beta0 <- as.numeric(Fcon%*%theta0 + xhat)[1:length(x)]

  fit.con <- mod$family$linkinv(model.matrix(mod)%*%beta0)
  dev.con <- sum(mod$family$dev.resids(mod$y, fit.con, wt = mod$prior.weights))

  # Use deviance-based F-statistic:

  zs <- c(zs, sqrt((dev.con - dev.mle)/(dev.mle/(n - p))))

  # Profile in the positive direction.

  i <- 2
  while (tail(zs, 1) < z.max){
    DzDpsi <- (zs[i] - zs[i - 1])/(ms[i] - ms[i-1])

    dpsi <- D.z / DzDpsi

    ms <- c(ms, tail(ms, 1) + dpsi)

    xhat <- qr.out$Q[, 1]*(ms[length(ms)]/qr.out$R[1])

    prof.out <- profile.lik(ms[length(ms)], theta0)

    pll <- prof.out$value
    theta0 <- prof.out$theta.con

    beta0 <- as.numeric(Fcon%*%theta0 + xhat)[1:length(x)]

    fit.con <- mod$family$linkinv(model.matrix(mod)%*%beta0)
    dev.con <- sum(mod$family$dev.resids(mod$y, fit.con, wt = mod$prior.weights)) ## Not sure about wt = 1 here.

    zs <- c(zs, sqrt((dev.con - dev.mle)/(dev.mle/(n - p))))

    i <- i + 1
  }

  # Profile in the negative direction.

  i <- 2
  while (zs[1] > -z.max){
    DzDpsi <- (zs[i-1] - zs[i])/(ms[i-1] - ms[i])

    dpsi <- D.z / DzDpsi

    ms <- c(ms[1] - dpsi, ms)

    xhat <- qr.out$Q[, 1]*(ms[1]/qr.out$R[1])

    prof.out <- profile.lik(ms[1], theta0)

    pll <- prof.out$value
    theta0 <- prof.out$theta.con

    beta0 <- as.numeric(Fcon%*%theta0 + xhat)[1:length(x)]

    fit.con <- mod$family$linkinv(model.matrix(mod)%*%beta0)
    dev.con <- sum(mod$family$dev.resids(mod$y, fit.con, wt = mod$prior.weights))

    zs <- c(-sqrt((dev.con - dev.mle)/(dev.mle/(n - p))), zs)

    i <- i + 1
  }

  resp.vals <- mod$family$linkinv(ms)

  profile <- approxfun(resp.vals, zs)

  # Not sure why we have to do this? Why aren't the resp.values always increasing?

  pconf <- if (resp.vals[2] - resp.vals[1] < 0) {
    function(m) 1-pt(profile(m), df = n - p)
  } else{
    function(m) pt(profile(m), df = n - p)
  }
  cconf <- function(m) pf(profile(m)^2, df1 = 1, df2 = n - p)

  dconf <-  function(m, dx = 1e-5) (pconf(m + dx) - pconf(m - dx))/(2*dx)

  qconf <- function(p) {
    fun.root <- function(z) pconf(z) - p

    rang.ms <- range(mod$family$linkinv(ms))

    return(uniroot(fun.root, interval = c(rang.ms[1], rang.ms[2]))$root) # Might need better upperbound here!
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(m) 1 - cconf(m)

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve)

  if (plot){
    plot.dconf(out, xlab = 'Expected Response')
    plot.cconf(out, conf.level = conf.level, xlab = 'Expected Response')
  }

  return(out)
}

#' Confidence Functions for the Coefficients of a GLM
#'
#' Confidence functions for the coefficients of a glm.
#'
#' @param mod an output from glm
#'
#' @return A list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each coefficient of the glm, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @examples
#' # Low birth weight example from page 38 of *Confidence, Likelihood, Probability*,
#'
#' data(lbw)
#'
#' mod <- glm(low ~ weight + age + black + other + smoker, data = lbw, family = binomial)
#'
#' mod.conf <- glm.beta.conf(mod)
#'
#' plot.cconf(mod.conf$weight)
#'
#' @export
glm.beta.conf <- function(mod){
  alpha = 1e-8
  prof <- profile(mod, alpha = alpha, maxsteps = 100, del = qnorm(1-alpha)/80)

  fam <- family(mod)

  Pnames <- names(B0 <- coef(mod))

  p <- length(Pnames)

  mf <- model.frame(mod)
  Y <- model.response(mf)
  n <- NROW(Y)

  switch(fam$family,
         binomial = ,
         poisson = ,
         `Negative Binomial` = {
           pivot.cdf <- function(x) pnorm(x)
           profName <- "z"
         }
         ,
         gaussian = ,
         quasi = ,
         inverse.gaussian = ,
         quasibinomial = ,
         quasipoisson = ,
         {
           pivot.cdf <- function(x) pt(x, n - p)
           profName <- "tau"
         }
  )

  return.list <- list()

  for (var.name in Pnames){
    return.list[[var.name]] <- list()

    z <- prof[[var.name]]$par.vals[, var.name]

    Hn <- pivot.cdf(prof[[var.name]][[profName]])

    # Confidence Distribution

    Cn <- approxfun(z, Hn)

    # Confidence Density

    dz <- min(diff(z))

    hn <- (Cn(z + dz) - Cn(z - dz))/(2*dz)

    cn <- approxfun(z, hn)

    # Confidence Quantile

    Qn <- approxfun(Hn, z)

    # Confidence Curve

    ccn <- approxfun(z, abs(2*Hn - 1))

    # P-curve

    Ps <- 2*apply(cbind(Hn, 1 - Hn), 1, min)

    Pn <- approxfun(z, Ps)

    # S-curve

    Sn <- approxfun(z, -log2(Ps))

    return.list[[var.name]][['pconf']] <- Cn
    return.list[[var.name]][['dconf']] <- cn
    return.list[[var.name]][['qconf']] <- Qn
    return.list[[var.name]][['cconf']] <- ccn
    return.list[[var.name]][['pcurve']] <- Pn
    return.list[[var.name]][['scurve']] <- Sn
  }

  class(return.list) <- 'glm.beta.conf'

  return(return.list)
}
