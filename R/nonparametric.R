#' @export sign.conf
sign.conf <- function(x, plot = TRUE, conf.level = 0.95){
  x.sort <- sort(x)
  
  x.min <- x.sort[1]
  x.max <- x.sort[length(x)]
  
  n <- length(x)
  
  # See page 66 of Hollander, Wolfe, and Chicken for P-value
  # calculation.
  
  # Note: this only gives the correct right-sided P-value,
  # will be too small for left-sided P-value, since should be
  # P(C >= c), not P(C > c).
  
  # An alternative is to use the mid P-value.
  pconf <- function(mu){
    c <- sum(x <= mu)
    
    # pval <- pbinom(c, n, 0.5) - 0.5*dbinom(c, n, 0.5)
    pval <- pbinom(c, n, 0.5)
    
    return(pval)
  }
  
  pconf <- Vectorize(pconf)
  
  # cconf <- function(mu) abs(2*pconf(mu) - 1)
  
  cconf <- function(mu) {
    c <- sum(x <= mu)
    
    if (c <= n/2){
      return(1 - 2*pbinom(c-1, n, 1/2))
    }else{
      return(1 - 2*pbinom(n - c, n, 1/2))
    }
    
    
  }
  
  cconf <- Vectorize(cconf)
  
  # Using procedure from page 80 of Hollander, Wolfe, and Chicken:
  
  qconf <- function(p) {
    # Using to make qconf give the standard two-sided interval for the sign test.
    if (p >= 0.5){
      x.sort[qbinom(1-p, size = n, prob = 1/2, lower.tail = FALSE) + 1]
    }else{
      x.sort[n + 1 - (qbinom(p, size = n, prob = 1/2, lower.tail = FALSE) + 1)]
    }
  }
  
  qconf <- Vectorize(qconf)
  
  pcurve <- function(mu) 1 - cconf(mu)
  
  scurve <- function(mu) -log2(pcurve(mu))
  
  out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)
  
  if (plot){
    plot.cconf(out, conf.level = conf.level, xlab = 'Median')
    
    cat(sprintf("Actual confidence level: %g\n\n", cconf(qconf((1 - conf.level)/2))))
  }
  
  return(out)
}

#' @export wilcox.conf
wilcox.conf <- function(x, y = NULL, plot = TRUE, conf.level = 0.95){
  if (is.null(y)){
    pconf <- Vectorize(function(mu) wilcox.test(x, mu = mu, alternative = 'greater')$p.value)
    
    cconf <- function(mu) abs(2*pconf(mu) - 1)
    
    qconf <- Vectorize(function(p) wilcox.test(x, alternative = 'less', conf.level = p, conf.int = TRUE)$conf.int[2])
    
    pcurve <- function(mu) 1 - cconf(mu)
    
    scurve <- function(mu) -log2(pcurve(mu))
    
    out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)
    
    if (plot){
      plot.cconf(out, conf.level = conf.level, xlab = '(Psuedo)Median')
    }
  }  else{
      pconf <- Vectorize(function(mu) wilcox.test(x, y, mu = mu, alternative = 'greater')$p.value)
      
      cconf <- function(mu) abs(2*pconf(mu) - 1)
      
      qconf <- Vectorize(function(p) wilcox.test(x, y, alternative = 'less', conf.level = p, conf.int = TRUE)$conf.int[2])
      
      pcurve <- function(mu) 1 - cconf(mu)
      
      scurve <- function(mu) -log2(pcurve(mu))
      
      out <- list(pconf = pconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)
      
      if (plot){
        plot.cconf(out, conf.level = conf.level, xlab = 'Shift')
      }
  }
  
  return(out)
}