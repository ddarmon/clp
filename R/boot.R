#' @export
bootcurve <- function(data, statistic, B = 2000, stratified = FALSE){
  if (stratified){
    strata <- data[, ncol(data)]

    boot.out <- boot(data = data, statistic = statistic, strata = strata, R = B)
  }else{
    boot.out <- boot(data = data, statistic = statistic, R = B)
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

confdist = function(bc, theta, param){
  Gn = bc$Gn[[param]]
  Phi.invs = qnorm(Gn(theta))

  # The BCa confidence distribution
  Hn = pnorm((Phi.invs - bc$z0[param])/(1 + bc$a[param]*(Phi.invs - bc$z0[param])) - bc$z0[param])

  return(Hn)
}

confdens = function(bc, param){
  # density.out = density(bc$t[, param], bw = "SJ") # Seems to undersmooth
  density.out = density(bc$t[, param], bw = "bcv", n = 1024) # Seems to oversmooth, which in this case is good.

  gn.percentile = density.out$y
  thetas = density.out$x

  z0 = bc$z0[param]; a = bc$a[param]

  Gn = bc$Gn[[param]]

  w = function(theta, Gn, z0, a){
    ztheta = qnorm(Gn(theta)) - z0

    bca.fac = dnorm(ztheta/(1+a*ztheta) - z0)/((1 + a*ztheta)^2*dnorm(ztheta + z0))

    return(bca.fac)
  }

  gn.bca = gn.percentile*w(thetas, Gn, z0, a)

  # Reweight so sums to 1, at the given discretization of theta.
  gn.bca = gn.bca/(sum(gn.bca, na.rm = TRUE)*diff(thetas[1:2]))

  dconf <- splinefun(thetas, gn.percentile)

  return(dconf)
}

confquant <- function(bc, p, param){
  Gn <- bc$Gn[[param]]
  z0 <- bc$z0[param]

  a <- bc$a[param]

  z <- qnorm(p)

  zp <- z0 + (z0 + z)/(1 - a*(z0 + z))

  Qn <- quantile(bc$t[, param], pnorm(zp))

  return(Qn)
}

t.two.sample <- function(data, id = 1:nrow(data)){
  dat <- data[id, ]

  d <- mean(dat[dat[, 2] == 1, 1]) - mean(dat[dat[, 2] == 2, 1])

  return(d)
}

t.one.sample <- function(data, id = 1:length(data)){
  dat <- data[id]

  d <- mean(dat)

  return(d)
}

#' @export t.boot.conf
t.boot.conf <- function(x, y = NULL, B = 2000, plot = TRUE, conf.level = 0.95){

  if (is.null(y)){
    dat <- x
    bc <- bootcurve(dat, t.one.sample, B = B)

    xlab <- 'mean'
  }else{
    dat <- cbind(c(x, y), c(rep(1, length(x)), rep(2, length(y))))
    bc <- bootcurve(dat, t.two.sample, B = B, stratified = TRUE)

    xlab <- 'mean[1] - mean[2]'
  }


  pconf <- function(x) confdist(bc, x, 1)
  cconf <- function(x) abs(2*pconf(x) - 1)

  # Only valid within the approximate range of the data!
  dconf <- function(x) confdens(bc, 1)(x)

  qconf <- function(p) confquant(bc, p, 1)

  pcurve <- function(x) 1 - cconf(x)
  scurve <- function(x) -log2(pcurve(x))

  out <- list(pconf = pconf, cconf = cconf, dconf = dconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    plot.dconf(out, xlab = xlab)
    plot.cconf(out, conf.level = conf.level, xlab = xlab)
  }

  return(out)
}