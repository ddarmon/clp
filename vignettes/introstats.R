## ---- message = FALSE---------------------------------------------------------
library(clp)

## -----------------------------------------------------------------------------
data(dietstudy)

head(dietstudy)

## ---- message = FALSE---------------------------------------------------------
library(ggformula)

gf_histogram(~ weightchange | diet, data = dietstudy)

## -----------------------------------------------------------------------------
conf.mu.lowcarb <- t.conf(dietstudy$weightchange[dietstudy$diet == 'Low Carb'])
conf.mu.lowfat <- t.conf(dietstudy$weightchange[dietstudy$diet == 'Low Fat'])

## -----------------------------------------------------------------------------
conf.diff <- t.conf(x = dietstudy$weightchange[dietstudy$diet == 'Low Fat'], y = dietstudy$weightchange[dietstudy$diet == 'Low Carb'])

## -----------------------------------------------------------------------------
conf.diff$qconf(c(0.025, 0.975))

## -----------------------------------------------------------------------------
conf.diff$pcurve(0)

