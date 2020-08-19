# Function to generate the multinomial confidence functions
# using bootstrapping and the method of Rudolf Beran (1988):
#
# Rudolf Beran. "Balanced simultaneous confidence sets." Journal of the American Statistical Association 83.403 (1988): 679-686.
make.multinomial.multicomp.obj <- function(x, n, K, Kinv){
  conf.out <- prop.conf(x, n, plot = FALSE)

  theta.hat <- x/n

  cconf <- function(theta) K(conf.out$cconf(theta))

  pconf <- function(theta){
    if (theta <= theta.hat){
      return(0.5 - 0.5*cconf(theta))
    }else{
      return(0.5 + 0.5*cconf(theta))
    }
  }

  pconf <- Vectorize(pconf)

  qconf <- function(p){
    fun.root <- function(x) pconf(x) - p

    uniroot(fun.root, interval = c(0, 1))$root
  }

  qconf <- Vectorize(qconf)

  pcurve <- function(theta) 1 - cconf(theta)

  scurve <- function(theta) -log2(pcurve(theta))

  out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  return(out)
}

#' Simultaneous Confidence Functions for Multinomial Proportions.
#'
#' Confidence functions for a multinomial proportions based on the
#' mid P-value for each proportion, with balancing across
#' proportions via bootstrapping with Beran's simultaneous
#' confidence sets.
#'
#' @param N a vector of counts in each category from a sample of size n
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#' @param B the number of bootstrap samples to use with Beran's method
#' @param col a vector of colors to use for plotting confidence functions of each proportion
#'
#' @return A list of lists of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each proportion without correction for multiple comparison (singlecomp) or with
#'         correction for multiple comparison (multicomp).
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, Likelihood, Probability. Vol. 41. Cambridge University Press, 2016.
#'
#'              Tore Schweder. "Confidence nets for curves." Advances In Statistical Modeling And Inference: Essays in Honor of Kjell A Doksum. 2007. 593-609.
#'
#'              Rudolf Beran. "Balanced simultaneous confidence sets." Journal of the American Statistical Association 83.403 (1988): 679-686.
#'
#' @examples
#' # A hypothetical Punnett square experiment with peas.
#'
#' N <- c(315, 108, 101, 32)
#'
#' names(N) <- c('Round, Yellow', 'Round, Green', 'Wrinkled, Yellow', 'Wrinkled, Green')
#'
#' col <- c('darkgoldenrod', 'darkgreen', 'darkgoldenrod1', 'darkolivegreen')
#'
#' peas.conf <- multinomial.conf(N, col = col)
#'
#' # Confidence intervals without and with correction for multiple comparisons:
#' peas.conf$singlecomp$`Round, Yellow`$qconf(c(0.025, 0.975))
#' peas.conf$multicomp$`Round, Yellow`$qconf(c(0.025, 0.975))
#'
#' @export multinomial.conf
multinomial.conf <- function(N, plot = TRUE, conf.level = 0.95, B = 2000, col = NULL){
  nam <- names(N)

  if (is.null(nam)){
    nam <- as.character(1:length(N))
  }

  n <- sum(N)

  alpha <- 1 - conf.level

  theta.mle <- N/sum(N)

  # (Number of bootstrap samples) x (Number of categories)
  boot.samps <- t(rmultinom(B, size = n, prob = theta.mle))

  V <- matrix(nrow = B, ncol = length(N))

  for (b in 1:B){
    cur.sample <- boot.samps[b, ]

    theta.boot <- cur.sample / sum(cur.sample)

    for (j in 1:length(N)){
      V[b, j] <- prop.conf(cur.sample[j], n, plot = FALSE)$cconf(theta.mle[j])
    }
  }

  Vmax <- apply(V, 1, max)
  K <- ecdf(Vmax)

  Kinv <- function(p) quantile(Vmax, p)

  singlecomp <- list()
  multicomp <- list()

  for (j in 1:length(N)){
    singlecomp[[nam[j]]] <- prop.conf(N[j], n, plot = FALSE)
    multicomp[[nam[j]]] <- make.multinomial.multicomp.obj(N[j], n, K, Kinv)
  }

  n.plot <- 2001

  if (plot){
    if (is.null(col)){
      col <- rainbow(length(N))
    }

    for (j in 1:length(N)){
      xlab <- sprintf('Proportion, Category %s', nam[j])

      plot.cconf(multicomp[[j]], xlab = xlab, col = col[j])
      curve(singlecomp[[j]]$cconf(x), lty = 2, add = TRUE, n = n.plot, col = col[j])
    }

    plot(0, 0, cex = 0, xlim = c(0, 1), ylim = c(0, 1.5), xlab = 'Proportions', ylab = 'Confidence Curve', yaxt = 'n')
    for (j in 1:length(N)){
      curve(multicomp[[j]]$cconf(x), add = TRUE, col = col[j], n = n.plot)
    }
    legend('topleft', legend = nam, lty = 1, col = col)
    axis(2, seq(0, 1, 0.25))
  }

  return(list(singlecomp = singlecomp, multicomp = multicomp))
}

#' Compute a P-value for Specified Multinomial Proportions
#'
#' Compute a P-value for specified multinomial proportions
#' using the confidence net returned by conf.multinomial().
#'
#' @param obj an output from conf.multinomial().
#' @param theta a vector as the same length as the N used
#'              to construct obj; the null value of the
#'              multinomial proportions.
#'
#' @return A P-value from theta.
#'
#' @examples
#' # A hypothetical Punnett square experiment with peas.
#'
#' N <- c(315, 108, 101, 32)
#'
#' names(N) <- c('Round, Yellow', 'Round, Green', 'Wrinkled, Yellow', 'Wrinkled, Green')
#'
#' col <- c('darkgoldenrod', 'darkgreen', 'darkgoldenrod1', 'darkolivegreen')
#'
#' peas.conf <- multinomial.conf(N, col = col)
#'
#' p.multinomial(peas.conf, c(9, 3, 3, 1)/16)
#'
#' @export p.multinomial
p.multinomial <- function(obj, theta){
  mod <- obj$multicomp

  pvals <- rep(NA, length(mod))

  for (j in 1:length(mod)){
    pvals[j] <- mod[[j]]$pcurve(theta[j])
  }

  return(min(pvals))
}

# <--------------------------------------------------------------->

ll.multinomial <- function(eta, gamma = NULL, N.flat){
  theta <- c(eta, 1 - sum(eta))

  # Negative, for minimization rather than maximization.
  return(-sum(N.flat*log(theta)))
}

grad.multinomial <- function(eta, gamma, N.flat){
  theta <- c(eta, 1 - sum(eta))

  # Additional term needed since each p[i, j] shows up
  # in the final parameter.

  minus.term <- N.flat[length(N.flat)]/theta[length(theta)]

  g <- N.flat / theta - minus.term

  # Negative, for minimization rather than maximization.
  return(-g[-length(g)])
}

# Compute chi-squared association parameter from J - 1 probabilities.
gamma.fun <- function(eta){
  theta <- c(eta, 1 - sum(eta))

  P <- matrix(theta, nrow = nrow(N))

  a <- rowSums(P)
  b <- colSums(P)

  # gamma.hat <- 0
  #
  # for (i in 1:nrow(N)){
  #   for (j in 1:ncol(N)){
  #     pr <- a[i]*b[j]
  #
  #     gamma.hat <- gamma.hat + (P[i, j] - pr)^2/pr
  #   }
  # }

  # Vectorized computation:

  ab <- outer(a, b)

  gamma.hat <- sum((P - ab)^2/ab)

  return(gamma.hat)
}

constraint.fun <- function(eta, gamma, N.flat){
  gamma.prof <- gamma.fun(eta)

  return(gamma.prof - gamma)
}

ineq.constraint.fun <- function(eta, gamma, N.flat){
  con <- c(eta, # ps > 0
           1 - eta) # ps < 1

  return(con)
}

ineq.constraint.jac <- function(eta, gamma, N.flat){
  jac <- rbind(pracma::eye(length(eta)), # for ps > 0
               -pracma::eye(length(eta))) # for ps < 1

  return(jac)
}

#' @export chisq.conf
chisq.conf <- function(N, plot = TRUE, conf.level = 0.95){
  # Flatten matrix for easier computes.
  N.flat <- as.numeric(N)

  phat.flat <- N.flat / sum(N.flat)

  gamma.hat <- gamma.fun(phat.flat[-length(phat.flat)])

  cconf <- function(gamma){
    # Runs faster **without** constraining 0 <= pi <= 1.
    #suppressWarnings(opt.out <- auglag(par = phat.flat[-length(phat.flat)], fn = ll.multinomial, gr = grad.multinomial, hin = ineq.constraint.fun, hin.jac = ineq.constraint.jac, heq = constraint.fun, control.outer = list(trace = FALSE, method = 'nlminb'), gamma = gamma, N.flat = N.flat))

    suppressWarnings(opt.out <- auglag(par = phat.flat[-length(phat.flat)], fn = ll.multinomial, gr = grad.multinomial, heq = constraint.fun, control.outer = list(trace = FALSE, method = 'nlminb'), gamma = gamma, N.flat = N.flat))

    return(pchisq(2*(ll.multinomial(opt.out$par, N.flat = N.flat) - ll.multinomial(phat.flat[-length(phat.flat)], N.flat = N.flat)), df = 1))
  }

  cconf <- Vectorize(cconf)

  pconf <- function(gamma){
    if (gamma < 0){
      return(0)
    }else{
      # Runs faster **without** constraining 0 <= pi <= 1.
      #suppressWarnings(opt.out <- auglag(par = phat.flat[-length(phat.flat)], fn = ll.multinomial, gr = grad.multinomial, hin = ineq.constraint.fun, hin.jac = ineq.constraint.jac, heq = constraint.fun, control.outer = list(trace = FALSE, method = 'nlminb'), gamma = gamma, N.flat = N.flat))

      suppressWarnings(opt.out <- auglag(par = phat.flat[-length(phat.flat)], fn = ll.multinomial, gr = grad.multinomial, heq = constraint.fun, control.outer = list(trace = FALSE, method = 'nlminb'), gamma = gamma, N.flat = N.flat))

      return(pnorm(sign(gamma - gamma.hat)*sqrt(2*(ll.multinomial(opt.out$par, N.flat = N.flat) - ll.multinomial(phat.flat[-length(phat.flat)], N.flat = N.flat)))))
    }
  }

  pconf <- Vectorize(pconf)

  qconf <- function(p){
    if (p <= pconf(0)){
      return(0)
    }else{
      fun.root <- function(x) pconf(x) - p

      right <- 0.1

      while (fun.root(right) < 0){
        right <- right*2
      }

      return(uniroot(fun.root, interval = c(0, right))$root)
    }
  }

  qconf <- Vectorize(qconf)

  dconf <- function(gamma, dx = 1e-5){
    if (gamma == 0){
      return(NA)
    }else{
      return((pconf(gamma + dx) - pconf(gamma - dx))/(2*dx))
    }
  }

  dconf <- Vectorize(dconf)

  pcurve <- function(gamma) 1 - cconf(gamma)

  scurve <- function(gamma) -log2(scurve(gamma))

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    xlab <- "Chi-square Parameter of Association"

    plot.dconf(out, xlab = xlab, n.points = 201)
    plot.cconf(out, xlab = xlab, n.points = 201, conf.level = conf.level)
  }

  return(out)
}
