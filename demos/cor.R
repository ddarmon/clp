source('../R/cor.R')
source('../R/plot.R')

load('../data/fat.rda')

fat <- fat[1:101, ]

# conf.out <- cor.conf(fat$body.fat, fat$weight, exact = FALSE)
conf.out <- cor.conf(fat$body.fat, fat$weight, exact = TRUE, plot = TRUE)

curve(conf.out$pconf(x), xlim = c(-1, 1), n = 36) # Works
curve(conf.out$pconf(x), xlim = c(-1, 1), n = 37) # Doesn't work

xs <- seq(-1, 1, length.out = 37)

# for (x in xs){
#   show(x)
#   conf.out$pconf(x)
# }

conf.out$pconf(xs[36])

dcorr <- function(x, rho, n){
  prefactor <- function(r, rho, n){
    (n - 2)*(1 - rho^2)^((n-1)/2)*(1 - r^2)^((n-4)/2)/pi
  }

  integrand <- function(w, r, rho, n){
    1/(cosh(w) - rho*r)^(n-1)
  }

  definite.integral <- function(r, rho, n){
    integrate(integrand, lower = 0, upper = Inf, r = r, rho = rho, n = n)$value
  }
  definite.integral <- Vectorize(definite.integral, vectorize.args = 'r')

  I <- definite.integral(x, rho, n)

  fn <- prefactor(x, rho, n)*I

  return(fn)
}

pcorr(cor(fat$body.fat, fat$weight), rho = xs[36], n = nrow(fat))
