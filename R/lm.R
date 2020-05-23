confcurve.lincom.lm <- function(mod, x){
  # Estimated expected response.
  mhat <- sum(x*coef(mod))

  # Estimated standard deviation of the
  # estimated expected response.
  s.mhat <- sqrt(as.numeric(t(x)%*%vcov(mod)%*%x))

  # Confidence distribution and curve based on
  # Gaussian distribution for errors.
  cd <- function(m) pt((m - mhat)/s.mhat, df = mod$df.residual)
  cc <- function(m) abs(2*cd(m) - 1)

  # Set up plotting limits.
  q.lim <- c(-1, 1)*qt(0.9999, df = mod$df.residual)
  m.lim <- mhat + q.lim*s.mhat

  # Generate plots.
  par(mfrow = c(1, 2))
  curve(cc, from = m.lim[1], to = m.lim[2], n = 2000,
        xlab = 'Expected Response', ylab = 'cc')
  abline(h = 0.95, lty = 3)
  curve(cd, from = m.lim[1], to = m.lim[2], n = 2000,
        xlab = 'Expected Response', ylab = 'cd')

  return(list(cc = cc, cd = cd))
}
