confcurve.lincom <- function(mod, x){
  # PDF / PMF from GLM
  dmod <- get_dmodel_function(mod)
  
  # Derivative of log-likelihood
  score <- get_score_function(mod) 
  
  # Log-likelihood evaluated at beta
  ell.beta <- function(beta){
    as.numeric(sum(log(dmod(coefficients = beta))))
  }
  
  score.beta <- function(beta){
    as.numeric(score(coefficients = beta))
  }
  
  beta.mle <- coef(mod)
  
  m.mle <- sum(x*beta.mle)
  
  # Get out feasible direction for the
  # constraint x*beta = m via QR decomp.
  
  qr.out <- householder(t(t(x)))
  
  Fcon <- qr.out$Q[, 2:length(x)]
  
  # Log-likelihood in terms of the feasible
  # direction.
  
  ell.con <- function(beta, xhat){
    beta <- matrix(beta, nrow = length(beta))
    
    arg <- Fcon%*%beta + xhat
    return(ell.beta(arg))
  }
  
  # Score in terms of the feasible direction.
  
  score.con <- function(beta, xhat){
    beta <- matrix(beta, nrow = length(beta))
    
    arg <- Fcon%*%beta + xhat
    
    t(Fcon)%*%score(coefficients = arg)
  }
  
  # The profile likelihood function.
  
  profile.lik <- function(m, beta0){
    opt.out <- optim(beta0, ell.con, score.con, method = 'BFGS', xhat = xhat, control = list(fnscale = -1))
    value <- opt.out$value
    beta.con <- opt.out$par
    
    return(list(value = value, beta.con = beta.con))
  }
  
  ms <- seq(m.mle - 2, m.mle + 2, length.out = 21)
  
  pll <- rep(0, length(ms))
  
  beta0 <- beta.mle[2:length(beta.mle)]
  
  for(m.ind in 1:length(ms)){
    m <- ms[m.ind]
    
    xhat <- qr.out$Q[, 1]*(m/qr.out$R[1])
    
    prof.out <- profile.lik(m, beta0)
    
    pll[m.ind] <- prof.out$value
    
    beta0 <- prof.out$beta.con
  }
  
  profile.lik.fun <- splinefun(mod$family$linkinv(ms), pll)
  
  dev.fun <- function(m){
    2*(ell.beta(beta.mle) - profile.lik.fun(m))
  }
  
  resp.vals <- mod$family$linkinv(ms)
  
  par(mfrow = c(1, 2))
  curve(profile.lik.fun, from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'Profile Likelihood')
  curve(dev.fun,  from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'Deviance')
  abline(h = qchisq(c(0.5, 0.9, 0.95, 0.99), 1), lty = 3)
  
  signed.sqrt.dev.fun <- function(m) sign(m - mod$family$linkinv(m.mle))*sqrt(dev.fun(m))
  
  cc <- function(m) pchisq(dev.fun(m), 1)
  cd <- function(m) pnorm(signed.sqrt.dev.fun(m))
  
  curve(cc(x), from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'cc')
  abline(h = 0.95, lty = 3)
  
  curve(cd(x), from = min(resp.vals), to = max(resp.vals),
       xlab = 'Expected Response', ylab = 'cd')
  
  return(list(cc = cc, cd = cd, profile.lik = profile.lik.fun, deviance = dev.fun))
}