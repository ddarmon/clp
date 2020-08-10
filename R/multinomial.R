make.multinomial.multicomp.fun <- function(x, n, K, Kinv){
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
  
  conf.singlecomp <- list()
  conf.multicomp <- list()
  
  for (j in 1:length(N)){
    conf.singlecomp[[j]] <- prop.conf(N[j], n, plot = FALSE)
    conf.multicomp[[j]] <- make.multinomial.multicomp.fun(N[j], n, K, Kinv)
  }
  
  n.plot <- 2001
  
  if (plot){
    if (is.null(col)){
      col <- rainbow(length(N))
    }
    
    for (j in 1:length(N)){
      xlab <- sprintf('Proportion, Category %s', nam[j])
      
      plot.cconf(conf.multicomp[[j]], xlab = xlab, col = col[j])
      curve(conf.singlecomp[[j]]$cconf(x), lty = 2, add = TRUE, n = n.plot, col = col[j])
    }
    
    plot(0, 0, cex = 0, xlim = c(0, 1), ylim = c(0, 1.5), xlab = 'Proportions', ylab = 'Confidence Curve', yaxt = 'n')
    for (j in 1:length(N)){
      curve(conf.multicomp[[j]]$cconf(x), add = TRUE, col = col[j], n = n.plot)
    }
    legend('topleft', legend = nam, lty = 1, col = col)
    axis(2, seq(0, 1, 0.25))
  }
  
  return(list(conf.singlecomp = conf.singlecomp, conf.multicomp = conf.multicomp))
}

#' @export p.multinomial
p.multinomial <- function(obj, theta){
  mod <- obj$conf.multicomp
  
  pvals <- rep(NA, length(mod))
  
  for (j in 1:length(mod)){
    pvals[j] <- mod[[j]]$pcurve(theta[j])
  }
  
  return(min(pvals))
}