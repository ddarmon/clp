cubic.root <- function(x, n, Delta){
  x0 <- x[1]; x1 <- x[2]
  n0 <- n[1]; n1 <- n[2]

  x <- x0 + x1
  n <- n0 + n1

  a0 <- x0*Delta*(1-Delta)
  a1 <- (n0*Delta - n - 2*x0)*Delta + x
  a2 <- (n1 + 2*n0)*Delta - n - x
  a3 <- n

  q <- a2^3/(3*a3)^3 - a1*a2/(6*a3^2) + a0/(2*a3)
  p <- sign(q)*sqrt(a2^2/(3*a3)^2 - a1 / (3*a3))
  a <- (1/3)*(pi + acos(q/p^3))

  p.con <- 2*p*cos(a) - a2/(3*a3)

  if (is.na(p.con)) { # Analytical solution sometimes fails, so resort to numerical solution.
    #   cat(sprintf("%g x^3 + %g x^2 + %g x + %g", a3, a2, a1, a0))
    roots <- polyroot(c(a0, a1, a2, a3))
    reroots <- Re(roots)

    if (Delta > 0){
      p.con <- reroots[reroots >= 0 & reroots <= 1 - Delta][1] # Select first, since sometimes get repeated roots.
    } else{
      p.con <- reroots[reroots >= -Delta & reroots <= 1][1]
    }
  }

  return(p.con)
}

cubic.root <- Vectorize(cubic.root, vectorize.args = c('Delta'))

#' @export
prop.conf <- function(x, n, plot = TRUE, conf.level = 0.95){
  if  (length(x) == 1 && length(n == 1)){
    out <- prop.conf.1s(x, n, plot = plot, conf.level = conf.level)
  } else if (length(x) == 2 && length(n) == 2){
    out <- prop.conf.2s(x, n, plot = plot, conf.level = conf.level)
  }

  return(out)
}

prop.conf.1s <- function(x, n, plot = TRUE, conf.level = 0.95){
  pconf <- function(p){
    cd <- vector(length = length(p))

    cd[p < 0] <- 0
    cd[p > 1] <- 1

    in.bounds <- (p >= 0) & (p <= 1)

    cd[in.bounds] <- pbinom(x, n, p[in.bounds], lower.tail = FALSE) + 0.5*dbinom(x, n, p[in.bounds])

    return(cd)
  }

  cconf <- function(p) abs(2*pconf(p) - 1)

  dconf <- function(p){
    if (x == 0){
      0.5*dbeta(p, 1, n)
    } else if (x == n){
      0.5*dbeta(p, n, 1)
    } else{
      0.5*dbeta(p, x, n-x+1) + 0.5*dbeta(p, x+1, n-x)
    }
  }

  qconf <- function(p){
    if (x == 0 && p <= 0.5){
      return(0)
    } else if (x == n && p >= 0.5){
      return(1)
    }else{
      fun.root <- function(z) {
        diff <- pconf(z) - p

        return(diff)
      }

      return(uniroot(fun.root, interval = c(0, 1))$root)
    }
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(p) 1 - cconf(p)

  out <- list(pconf = pconf, dconf = dconf, qconf = qconf, cconf = cconf, pcurve = pcurve)

  if (plot){
    plot.dconf(out, xlab = 'Proportion p')
    plot.cconf(out, conf.level = conf.level, xlab = 'Proportion p')
  }

  return(out)
}

prop.conf.2s <- function(x, n, plot = TRUE, conf.level = 0.95){
  x0 <- x[1]; x1 <- x[2]
  n0 <- n[1]; n1 <- n[2]

  pconf.score <- function(Delta){
    if (Delta < -1){
      return(0)
    }else if (Delta == -1){
      Delta = -1+1e-5 # To handle potential atom at Delta = -1.
    }else if (Delta >= 1){
      return(1)
    }

    p0 <- x0/n0
    p1 <- x1/n1

    NUM <- Delta - (p0 - p1)

    p1.Delta <- cubic.root(c(x0, x1), c(n0, n1), Delta)
    p0.Delta <- p1.Delta + Delta

    D0 <- p0.Delta*(1-p0.Delta)/n0
    D1 <- p1.Delta*(1-p1.Delta)/n1

    fac <- (n0 + n1)/(n0 + n1 - 1) # Miettinen and Nurminen correction.

    DENOM <- sqrt(fac*(D0 + D1))

    return(pnorm(NUM/DENOM))
  }

  pconf.score <- Vectorize(pconf.score)

  cconf.score <- function(Delta) abs(2*pconf.score(Delta) - 1)

  dconf.score <- function(Delta, dx = 1e-10){
    dd <- (pconf.score(Delta + dx) - pconf.score(Delta - dx))/(2*dx)

    dd[(Delta == -1) | (Delta == 1)] <- 0

    return(dd)
  }

  qconf.score <- function(p){
    if ((x0 == 0) && (x1 == n1) && (p <= 0.5)){ # Handle atom at Delta = -1.
      return(-1)
    }else{
      fun.root <- function(z) {
        diff <- pconf.score(z) - p

        return(diff)
      }

      return(uniroot(fun.root, interval = c(-1, 1))$root)
    }
  }

  qconf.score <- Vectorize(qconf.score)

  pcurve.score <- function(Delta) 1 - cconf.score(Delta)

  out <- list(pconf = pconf.score, dconf = dconf.score, qconf = qconf.score, cconf = cconf.score, pcurve = pcurve.score)

  if (plot){
    plot.dconf(out, xlab = 'Difference (p[1] - p[2])')
    plot.cconf(out, conf.level = conf.level, xlab = 'Difference (p[1] - p[2])')
  }

  return(out)
}

prop.conf.agresti.caffo <- function(x, n, plot = TRUE, conf.level = 0.95){
  x0 <- x[1]; x1 <- x[2]
  n0 <- n[1]; n1 <- n[2]

  x0 <- x0 + 1
  x1 <- x1 + 1

  n0 <- n0 + 2
  n1 <- n1 + 2

  p0 <- x0/n0
  p1 <- x1/n1

  pconf.ac <- function(Delta){
    NUM <- Delta - (p0 - p1)

    DENOM <- sqrt(p0*(1-p0)/n0 + p1*(1-p1)/n1)

    return(pnorm(NUM/DENOM))
  }

  dconf.ac <- function(Delta){
    NUM <- Delta - (p0 - p1)

    DENOM <- sqrt(p0*(1-p0)/n0 + p1*(1-p1)/n1)

    return(dnorm(NUM/DENOM)/DENOM)
  }

  cconf.ac <- function(Delta) abs(2*pconf.ac(Delta) - 1)

  qconf.ac <- function(p){
    DENOM <- sqrt(p0*(1-p0)/n0 + p1*(1-p1)/n1)

    return((p0 - p1) + qnorm(p)*DENOM)
  }

  pcurve.ac <- function(Delta) 1 - cconf.ac(Delta)

  out <- list(pconf = pconf.ac, dconf = dconf.ac, qconf = qconf.ac, cconf = cconf.ac, pcurve = pcurve.ac)

  if (plot){
    plot.dconf(out, xlab = 'Difference (p[1] - p[2])')
    plot.cconf(out, conf.level = conf.level, xlab = 'Difference (p[1] - p[2])')
  }

  return(out)
}

#' @export
risk.conf <- function(x, n, plot = TRUE, conf.level = 0.95, log = ''){
  x0 <- x[1]; x1 <- x[2]
  n0 <- n[1]; n1 <- n[2]



  if ((x0 == 0) && (x1 == 0)){
    warning("Confidence functions not defined when no successes in either group.")
    return(FALSE)
  }

  pconf.score <- function(rho){
    if (rho <= 0){
      return(0)
    }

    p0 <- x0/n0
    p1 <- x1/n1

    x <- (x0 + x1)
    n <- (n0 + n1)

    NUM <- (p1 - p0*rho)^2

    A <- n*rho
    B <- -(n1*rho + x1 + n0 + x0*rho)
    C <- x

    p0.rho <- (-B - sqrt(B^2 - 4*A*C))/(2*A)
    p1.rho <- p0.rho*rho

    S0 <- rho^2*p0.rho*(1-p0.rho)/n0
    S1 <- p1.rho*(1-p1.rho)/n1

    fac <- n/(n-1)

    DENOM <- (S0 + S1)*fac

    return(pnorm(sign(p0*rho - p1)*sqrt(NUM/DENOM)))
  }

  pconf.score <- Vectorize(pconf.score)

  cconf.score <- function(rho) abs(2*pconf.score(rho) - 1)

  dconf.score <- function(rho, dx = 1e-5){
    dd <- (pconf.score(rho + dx) - pconf.score(rho - dx))/(2*dx)

    if ((x1 == 0) || (x0 == n0)){
      dd[rho == 0] = 0
    }

    return(dd)
  }

  qconf.score <- function(p){
    # Need to deal with atoms:

    # When x1 = 0, atom w/ prob 1/2 at rho = 0

    # When x0 = 0, atom w/ prob 1/2 at rho = Infty

    if (x1 == 0 && p <= 0.5){
      return(0)
    }else if (x0 == 0 && p >= 0.5){
      return(Inf)
    }

    fun.root <- function(z) {
      rr <- pconf.score(z) - p

      return(rr)
    }

    ub <- 100

    while (pconf.score(ub) < p){
      ub <- ub*10
    }

    return(uniroot(fun.root, interval = c(0, ub))$root) # Might need better upperbound here!
  }

  qconf.score <- Vectorize(qconf.score)

  pcurve.score <- function(rho) 1 - cconf.score(rho)

  out <- list(pconf = pconf.score, dconf = dconf.score, qconf = qconf.score, cconf = cconf.score, pcurve = pcurve.score)

  if (plot){
    xlim <- qconf.score(c(0.001, 0.999))

    # For Relative risk atoms:
    if (xlim[2] == Inf){
      xlim[2] <- qconf.score(0.4)
    }

    if (log == 'x'){
      xlab.cconf = 'Relative Risk (p[2]/p[1]), Log-scaling'
    }else{
      xlab.cconf = 'Relative Risk (p[2]/p[1])'
    }

    plot.dconf(out, xlab = 'Relative Risk (p[2]/p[1])', xlim = xlim)
    plot.cconf(out, conf.level = conf.level, xlab = xlab.cconf, xlim = xlim, log = log)
  }

  return(out)
}

#' @export
oddsratio.conf <- function(x, n, plot = TRUE, conf.level = 0.95, log = ''){
  ## See page 236 of *Confidence, Likelihood, Probability*
  ## for the form of the confidence distribution
  ## for the odds ratio.

  y1 <- x[1]; y2 <- x[2]
  n1 <- n[1]; n2 <- n[2]

  p1 <- y1/n1
  p2 <- y2/n2

  or <- (p2*(1-p1))/(p1*(1-p2))

  # Let rho be the odds ratio, and
  #     psi be the log-odds ratio:
  #
  # psi = log(rho)

  log.or <- log(or)

  kappa.s <- 1/y1 + 1/(n1 - y1) + 1/y2 + 1/(n2 - y2)

  m1 <- y1 + y2

  # Use normal approximation to determine a
  # reasonable upper-bound for the odds ratio
  # rho:

  log.or.upper <- qnorm(0.99999, mean = log.or, sd = sqrt(kappa.s))

  or.upper <- exp(log.or.upper)

  n <- 1000
  rhos <- seq(0.01, or.upper, length.out = n)

  pconf <- function(rho){
    if (rho <= 0){
      return(0)
    }else{
      p.i <- dnoncenhypergeom(x = NA, n1 = n1, n2 = n2, m1 = m1, psi = rho)

      sum.inds <- which(p.i[, 1] > y2)

      p.i <- p.i[, 2]

      C <- sum(p.i[sum.inds]) + 0.5*p.i[sum.inds[1] - 1]

      return(C)
    }
  }

  pconf <- Vectorize(pconf)

  cconf <- function(rho) abs(2*pconf(rho) - 1)
  pcurve <- function(rho) 1 - cconf(rho)
  scurve <- function(rho) -log2(pcurve(rho))

  dconf <- function(rho, dx = 1e-10) (pconf(rho + dx) - pconf(rho - dx))/(2*dx)

  qconf <- function(p){
    fun.root <- function(z) {
      diff <- pconf(z) - p

      return(diff)
    }

    while (fun.root(or.upper) < 0){
      or.upper <- or.upper*10
    }

    return(uniroot(fun.root, interval = c(0, or.upper))$root)
  }

  qconf <- Vectorize(qconf)

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    xlim <- qconf(c(0.001, 0.999))

    if (log == 'x'){
      xlab.cconf = 'Odds Ratio (odds[2]/odds[1]), Log-scaling'
    }else{
      xlab.cconf = 'Odds Ratio (odds[2]/odds[1])'
    }

    plot.dconf(out, xlab = 'Odds Ratio (odds[2]/odds[1])', xlim = xlim)
    plot.cconf(out, conf.level = conf.level, xlab = xlab.cconf, xlim = xlim, log = log)
  }

  return(out)
}
