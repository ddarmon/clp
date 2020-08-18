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
#' plot.cconf(beta.conf$weight)
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
