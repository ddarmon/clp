#' @export
lm.lincom.conf <- function(mod, x, plot = TRUE, conf.level = 0.95){
  # Estimated expected response.
  mhat <- sum(x*coef(mod))

  # Estimated standard deviation of the
  # estimated expected response.
  s.mhat <- sqrt(as.numeric(t(x)%*%vcov(mod)%*%x))

  # Confidence function based on
  # Gaussian distribution for errors.
  pconf <- function(m) pt((m - mhat)/s.mhat, df = mod$df.residual)
  dconf <- function(m) dt((m - mhat)/s.mhat, df = mod$df.residual)/s.mhat
  cconf <- function(m) abs(2*pconf(m) - 1)
  qconf <- function(p) mhat + s.mhat*qt(p, df = mod$df.residual)

  pcurve <- function(m) 1 - cconf(m)

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve)

  if (plot){
    plot.dconf(out, xlab = 'Expected Response')
    plot.cconf(out, conf.level = conf.level, xlab = 'Expected Response')
  }

  return(out)
}

#' @export
lm.beta.conf <- function(mod){
  b <- coef(mod)
  s.b <- sqrt(diag(vcov(mod)))

  pconf <- function(var.name) return(function(beta) pt((beta - b[var.name])/s.b[var.name], df = mod$df.residual))
  dconf <- function(var.name) return(function(beta) dt((beta - b[var.name])/s.b[var.name], df = mod$df.residual)/s.b[var.name])
  cconf <- function(var.name) return(function(beta) abs(2*pconf(var.name)(beta)-1))

  pcurve <- function(var.name) return(function(beta) 1 - cconf(var.name)(beta))
  scurve <- function(var.name) return(function(beta) -log2(pcurve(var.name)(beta)))

  qconf <- function(var.name) return(function(p) b[var.name] + s.b[var.name]*qt(p, df = mod$df.residual))

  out <- list(variable.names = names(coef(mod)), pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  class(out) <- 'lm.beta.conf'

  return(out)
}
