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
  Pnames <- names(B0 <- coef(mod))

  p <- length(Pnames)

  S.b <- diag(vcov(mod.lm))

  return.list <- list()

  for (var.name in Pnames){
    return.list[[var.name]] <- list()

    b   <- B0[var.name]
    s.b <- sqrt(S.b[var.name])

    pconf <- function(beta) pt((beta - b)/s.b, df = mod$df.residual)
    dconf <- function(beta) dt((beta - b)/s.b, df = mod$df.residual)/s.b

    cconf <- function(beta) abs(2*pconf(beta) - 1)
    pcurve <- function(beta) 1-cconf(beta)
    scurve <- function(beta) -log2(pcurve(beta))


    qconf <- function(p) b + s.b*qt(p, df = mod$df.residual)

    return.list[[var.name]][['pconf']] <- pconf
    return.list[[var.name]][['dconf']] <- dconf
    return.list[[var.name]][['qconf']] <- qconf
    return.list[[var.name]][['cconf']] <- cconf
    return.list[[var.name]][['pcurve']] <- pcurve
    return.list[[var.name]][['scurve']] <- scurve
  }

  return(return.list)
}
