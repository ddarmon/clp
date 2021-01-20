#' Confidence Functions for the Expected Response of a Linear Model
#'
#' Confidence functions for the expected response of a linear model under
#' Gaussian noise assumptions.
#'
#' @param mod an output from lm
#' @param x a vector containing the predictor values for which the expected
#'          response is desired
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the expected response at x, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @examples
#' # Prediction of body fat percentage from other body measurements.
#'
#' data(fat)
#'
#' mod <- lm(body.fat ~ age + weight + height, data = fat)
#'
#' # Confidence curve for the expected BFP of a 170lb, 6ft, 20-year old male
#'
#' x <- c(1,    # Intercept
#'        20,   # 20 years old
#'        170,  # 170 lbs
#'        6*12) # 72 inches (6 feet)
#'
#' m.conf <- lm.lincom.conf(mod, x)
#'
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
    display.dconf(out, xlab = 'Expected Response')
    display.cconf(out, conf.level = conf.level, xlab = 'Expected Response')
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

#' Confidence Functions for the Coefficients of a Linear Model
#'
#' Confidence functions for the coefficients of a linear model under
#' Gaussian noise assumptions.
#'
#' @param mod an output from lm
#'
#' @return A list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each coefficient of the linear model, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @examples
#' # Prediction of body fat percentage from other body measurements.
#'
#' data(fat)
#'
#' mod <- lm(body.fat ~ age + weight + height, data = fat)
#'
#' beta.conf <- lm.beta.conf(mod)
#'
#' display.cconf(beta.conf$weight)
#'
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

  out <- betas

  class(out) <- 'lm.beta.conf'

  return(out)
}

#' Confidence Functions for the Noise Standard Deviation of a Linear Model
#'
#' Confidence functions for the noise standard deviation of a linear model under
#' Gaussian noise assumptions.
#'
#' @param mod an output from lm
#' @param plot whether to plot the confidence density and curve
#' @param conf.level the confidence level for the confidence interval indicated on the confidence curve
#'
#' @return A list containing the confidence functions pconf, dconf, cconf, and qconf
#'         for the noise standard deviation, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#' @examples
#' # Prediction of body fat percentage from other body measurements.
#'
#' data(fat)
#'
#' mod <- lm(body.fat ~ age + weight + height, data = fat)
#'
#' sigma.conf <- lm.sigma.conf(mod)
#'
#' display.cconf(sigma.conf)
#'
#' @export lm.sigma.conf
lm.sigma.conf <- function(mod, plot = TRUE, conf.level = 0.95){
  df <- mod$df.residual

  sigma.hat <- summary(mod)$sigma

  var.hat <- sigma.hat^2

  pconf <- function(sigma) pchisq(df*var.hat/sigma^2, df = df, lower.tail = FALSE)

  dconf <- function(sigma) dchisq(df*var.hat/sigma^2, df = df)*(2*df*var.hat)/(sigma^3)

  qconf <- function(p) sqrt(df*var.hat/qchisq(1-p, df))

  cconf <- function(sigma) abs(2*pconf(sigma) - 1)

  pcurve <- function(sigma) 1 - cconf(sigma)

  scurve <- function(sigma) -log2(pcurve(sigma))

  out <- list(pconf = pconf, dconf = dconf, cconf = cconf, qconf = qconf, pcurve = pcurve, scurve = scurve)

  if (plot){
    display.dconf(out, xlab = 'Noise Standard Deviation')
    display.cconf(out, conf.level = conf.level, xlab = 'Noise Standard Deviation')
  }

  return(out)
}
