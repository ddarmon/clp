confcurve.lincom <- function(mod, x){
  # PDF / PMF from GLM
  dmod <- get_dmodel_function(mod)

  # Derivative of log-likelihood
  score <- get_score_function(mod)

  # Log-likelihood evaluated at beta
  ell.beta <- function(beta){
    as.numeric(sum(log(dmod(coefficients = beta))))
  }

  beta.mle <- coef(mod)
  ell.mle <- ell.beta(beta.mle)

  m.mle <- sum(x*beta.mle)

  # Get out feasible direction for the
  # constraint x*beta = m via QR decomp.

  qr.out <- householder(t(t(x)))

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

  #Use procedure as described on page 23 of
  #
  # *Fitting Linear Mixed-Effects Models Using lme4*
  # (Bates, Machler, Bolker, and Walker)
  #
  # to perform the profiling.

  ms <- c(m.mle)
  zs <- c(0)

  z.max <- qnorm(0.999)
  D.z <- z.max / 16

  # Probably should take advantage of asymptotic results here
  # rather than use a fudge factor.

  new.m <- if (ms == 0) 0.001 else 1.01*ms

  ms <- c(ms, new.m)

  # Get value at initial new point:

  xhat <- qr.out$Q[, 1]*(ms[2]/qr.out$R[1])
  beta0 <- beta.mle[2:length(beta.mle)]
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

  par(mfrow = c(1, 2))
  curve(profile.lik.fun, from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'Profile Likelihood')
  curve(dev.fun,  from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'Deviance')
  abline(h = qchisq(c(0.5, 0.9, 0.95, 0.99), 1), lty = 3)

  signed.sqrt.dev.fun <- function(m) sign(m - mod$family$linkinv(m.mle))*sqrt(dev.fun(m))

  cc <- function(m) pchisq(dev.fun(m), 1)
  cd <- function(m) pnorm(signed.sqrt.dev.fun(m))

  curve(cc(x), from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'cc', n = 2001)
  abline(h = 0.95, lty = 3)

  curve(cd(x), from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'cd', n = 2001)

  return(list(cc = cc, cd = cd, profile.lik = profile.lik.fun, deviance = dev.fun))
}

confcurve.lincom.disp <- function(mod, x){
  n <- nrow(mod$model)
  p <- length(mod$coefficients)

  # PDF / PMF from GLM
  dmod <- get_dmodel_function(mod)

  # Derivative of log-likelihood
  score <- get_score_function(mod)

  # Log-likelihood evaluated at theta
  ell.theta <- function(theta){
    as.numeric(sum(log(dmod(coefficients = theta[1:(length(theta)-1)], dispersion = theta[length(theta)]))))
  }

  score.theta <- function(theta){
    as.numeric(score(coefficients = theta[1:(length(theta)-1)], dispersion = theta[length(theta)]))
  }

  # MLE of Dispersion Parameter
  phi.hat <- enrich(mod, with = "mle of dispersion")$dispersion_mle

  theta.mle <- c(coef(mod), phi.hat)

  ell.mle <- ell.theta(theta.mle)

  m.mle <- sum(x*theta.mle[1:length(x)])

  # Get out feasible direction for the
  # constraint x*theta = m via QR decomp.

  x.aug <- c(x, 0)

  qr.out <- householder(t(t(x.aug)))

  Fcon <- qr.out$Q[, 2:length(x.aug)]

  # Log-likelihood in terms of the feasible
  # direction.

  ell.con <- function(theta, xhat){
    theta <- matrix(theta, nrow = length(theta))
    theta[length(theta)] <- theta[length(theta)]^2

    arg <- Fcon%*%theta + xhat
    return(ell.theta(arg))
  }

  # Score in terms of the feasible direction.

  score.con <- function(theta, xhat){
    theta <- matrix(theta, nrow = length(theta))
    theta[length(theta)] <- theta[length(theta)]^2

    arg <- Fcon%*%theta + xhat

    t(Fcon)%*%score(coefficients = arg[1:(length(arg)-1)], dispersion = arg[length(arg)])
  }

  # The profile likelihood function.

  profile.lik <- function(m, theta0){
    opt.out <- optim(theta0, ell.con, score.con, method = 'BFGS', xhat = xhat, control = list(fnscale = -1))
    value <- opt.out$value
    theta.con <- opt.out$par

    return(list(value = value, theta.con = theta.con))
  }

  #Use procedure as described on page 23 of
  #
  # *Fitting Linear Mixed-Effects Models Using lme4*
  # (Bates, Machler, Bolker, and Walker)
  #
  # to perform the profiling.

  ms <- c(m.mle)
  zs <- c(0)

  z.max <- qt(0.999, n - p)
  D.z <- z.max / 16

  # Probably should take advantage of asymptotic results here
  # rather than use a fudge factor.

  new.m <- if (ms == 0) 0.001 else 1.01*ms

  ms <- c(ms, new.m)

  # Get value at initial new point:

  xhat <- qr.out$Q[, 1]*(ms[2]/qr.out$R[1])
  theta0 <- theta.mle[2:length(theta.mle)]
  prof.out <- profile.lik(ms[2], theta0)

  pll <- prof.out$value
  theta0 <- prof.out$theta.con

  zs <- c(zs, sqrt(2*(ell.mle - pll)))

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

    zs <- c(zs, sqrt(2*(ell.mle - pll)))

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

    zs <- c(-sqrt(2*(ell.mle - pll)), zs)

    i <- i + 1
  }

  # Recover profile log-likelihood from Z-values.

  pll <- ell.mle - (zs^2)/2

  profile.lik.fun <- splinefun(mod$family$linkinv(ms), pll)

  dev.fun <- function(m){
    2*(ell.mle - profile.lik.fun(m))
  }

  resp.vals <- mod$family$linkinv(ms)

  par(mfrow = c(1, 2))
  curve(profile.lik.fun, from = min(resp.vals), to = max(resp.vals),
        xlab = 'Expected Response', ylab = 'Profile Likelihood')
  curve(dev.fun,  from = min(resp.vals), to = max(resp.vals),
        xlab = 'Expected Response', ylab = 'Deviance')
  abline(h = qf(c(0.5, 0.9, 0.95, 0.99), df1 = 1, df2 = n - p), lty = 3)

  signed.sqrt.dev.fun <- function(m) sign(m - mod$family$linkinv(m.mle))*sqrt(dev.fun(m))

  cc <- function(m) pf(dev.fun(m), df1 = 1, df2 = n - p)
  cd <- function(m) pt(signed.sqrt.dev.fun(m), df = n - p)

  curve(cc(x), from = min(resp.vals), to = max(resp.vals),
        xlab = 'Expected Response', ylab = 'cc', n = 2001)
  abline(h = 0.95, lty = 3)

  curve(cd(x), from = min(resp.vals), to = max(resp.vals),
        xlab = 'Expected Response', ylab = 'cd', n = 2001)

  return(list(cc = cc, cd = cd, profile.lik = profile.lik.fun, deviance = dev.fun))
}
