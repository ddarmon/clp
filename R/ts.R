#' Confidence Functions for the Coefficients of an ARIMA Model
#'
#' Approximate confidence functions for the coefficients of an ARIMA model under
#' Gaussian noise assumptions.
#'
#' @param mod an output from stats::arima or astsa::sarima
#'
#' @return A list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each coefficient of the ARIMA model, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#' Peter J. Brockwell and Richard A. Davis. Introduction to Time Series and Forecasting. Springer, New York, 1996. Sections 3.3 and 8.3.
#'
#'
#' @examples
#' # Prediction of luteinizing hormone time series.
#'
#' data(lh)
#'
#' mod <- arima(lh, order = c(3,0,0))
#'
#' beta.conf <- arima.conf(mod)
#'
#' display.cconf(beta.conf$ar1)
#'
#' @export
arima.conf <- function(mod){
  if ("degrees_of_freedom" %in% names(mod)){ # mod is an output from astsa::sarima
    vnames <- rownames(mod$ttable)

    b <- mod$ttable[, 1]
    s.b <- mod$ttable[, 2]

    df <- mod$degrees_of_freedom
  }else{ # Assume mod is an output from stats::arima
    b <- mod$coef

    vnames <- names(b)

    s.b <- sqrt(diag(mod$var.coef))

    df <- mod$nobs-length(b)
  }

  betas <- list()

  for (var in vnames){
    betas[[var]] <- make.lm.beta.conf.obj(as.numeric(b[var]), as.numeric(s.b[var]), df)
  }

  out <- betas

  class(out) <- 'arima.conf'

  return(out)
}
