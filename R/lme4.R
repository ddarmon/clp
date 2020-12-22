#' Confidence Functions for the Coefficients of a LMM / GLMM
#'
#' Confidence functions for the Coefficients of a (Generalized) Linear Mixed Model.
#'
#' @param mod an output from glmer
#'
#' @return A list of lists containing the confidence functions pconf, dconf, cconf, and qconf
#'         for each coefficient of the GLMM, as well as the P-curve and S-curve.
#'
#' @references  Tore Schweder and Nils Lid Hjort. Confidence, likelihood, probability. Vol. 41. Cambridge University Press, 2016.
#'
#'
#' @examples
#'
#' @export
glmer.beta.conf <- function(mod){
  alpha = 1e-8
  # Producs a more finely-explored profile confidence distribution
  # prof <- profile(mod, alphamax = alpha, maxpts = 100, delta = qnorm(1-alpha)/80)
  prof <- profile(mod)

  fam <- family(mod)

  Pnames <- names(B0 <- fixef(mod))

  p <- length(Pnames)

  mf <- model.frame(mod)
  Y <- model.response(mf)
  n <- NROW(Y)

  switch(fam$family,
         binomial = ,
         poisson = ,
         `Negative Binomial` = {
           pivot.cdf <- function(x) pnorm(x)
           profName <- '.zeta'
         }
         ,
         gaussian = ,
         quasi = ,
         inverse.gaussian = ,
         quasibinomial = ,
         quasipoisson = ,
         {
           pivot.cdf <- function(x) pt(x, n - p)
           profName <- '.zeta'
         }
  )

  return.list <- list()

  for (var.name in Pnames){
    cur.inds <- prof$.par == var.name

    return.list[[var.name]] <- list()

    z <- prof[cur.inds, var.name]

    Hn <- pivot.cdf(prof[cur.inds, profName])

    # Confidence Distribution

    Cn <- approxfun(z, Hn)

    # Confidence Density

    dz <- min(diff(z))

    hn <- (Cn(z + dz) - Cn(z - dz))/(2*dz)

    cn <- approxfun(z, hn)

    # Confidence Quantile

    Qn <- approxfun(Hn, z)

    # Confidence Curve

    ccn <- approxfun(z, abs(2*Hn - 1))

    # P-curve

    Ps <- 2*apply(cbind(Hn, 1 - Hn), 1, min)

    Pn <- approxfun(z, Ps)

    # S-curve

    Sn <- approxfun(z, -log2(Ps))

    return.list[[var.name]][['pconf']] <- Cn
    return.list[[var.name]][['dconf']] <- cn
    return.list[[var.name]][['qconf']] <- Qn
    return.list[[var.name]][['cconf']] <- ccn
    return.list[[var.name]][['pcurve']] <- Pn
    return.list[[var.name]][['scurve']] <- Sn
  }

  class(return.list) <- 'glmer.beta.conf'

  return(return.list)
}
