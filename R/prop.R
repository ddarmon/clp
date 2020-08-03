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

prop.conf <- function(x, n, plot = TRUE){
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

    DENOM <- sqrt(D0 + D1)

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

  if (plot){
    tail.prob <- 0.001
    lims <- qconf.score(c(tail.prob/2, 1 - tail.prob/2))

    curve(dconf.score(x), xlim = lims, xlab = 'Difference (p[1] - p[2])', ylab = 'Confidence Density', n = 2001)
    curve(cconf.score(x), xlim = lims, xlab = 'Difference (p[1] - p[2])', ylab = 'Confidence Density', n = 2001)
    segments(x0 = qconf.score(0.025), x1 = qconf.score(0.975), y0 = 0.95, col = 'red')
  }

  out <- list(pconf = pconf.score, dconf = dconf.score, qconf = qconf.score, cconf = cconf.score)

  return(out)
}

prop.conf.agresti.caffo <- function(x, n){
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

  out <- list(pconf = pconf.ac, dconf = dconf.ac, qconf = qconf.ac, cconf = cconf.ac)
}
