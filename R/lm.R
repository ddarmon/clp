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

make.lm.beta.conf.obj <- function(b, s.b, df){
  b; s.b; df # Need to do this to instantiate these, so they save
             # after encapsulation.
  pconf <- function(beta) pt((beta - b)/s.b, df = df)
  dconf <- function(beta) dt((beta - b)/s.b, df = df)/s.b
  cconf <- function(beta) abs(2*pconf(beta)-1)

  pcurve <- function(beta) 1 - cconf(beta)
  scurve <- function(beta) -log2(pcurve(beta))

  qconf <- function(p) b + s.b*qt(p, df = df)

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  return(out)
}

#' @export
lm.beta.conf <- function(mod){
  vnames <- names(coef(mod))

  b <- coef(mod)
  s.b <- sqrt(diag(vcov(mod)))

  df <- mod$df.residual

  betas <- list()

  for (var in vnames){
    betas[[var]] <- make.lm.beta.conf.obj(as.numeric(b[var]), as.numeric(s.b[var]), df)
  }

  out <- list(variable.names = vnames, betas = betas)

  class(out) <- 'lm.beta.conf'

  return(out)
}
