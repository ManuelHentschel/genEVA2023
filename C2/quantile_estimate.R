################################################################
################################################################
################################################################
### GPD quantile estimate
## Auxiliary functions for GPD quatile estimation
################################################################
### 
library(extRemes)
quantile_estimate <- function(th1, margin, plot=F,return_fit=F,period=200){
  u   <- quantile(margin,th1)
  fit <- mev::fit.gpd(
    margin,
    threshold = u,
    method = "Grimshaw",
    show = FALSE,
    MCMC = NULL,
    k = 4,
    tol = 1e-08,
    fpar = NULL,
    warnSE = FALSE
  )
  if(plot) plot(fit)
  q0 <- evd::qgpd(1-1/(300*period*(1-th1)), 
                  scale = fit$estimate[1], 
                  shape =  fit$estimate[2] ) + u[1]
  q0 <- as.double(q0)
  if(return_fit) return(list(fit=fit,q0=q0) )
  else return(q0)
}
##
quantile_ci <- function(th1, margin, plot=F,return_fit=F,period=200,estimates=c(140,250)){
  u   <- quantile(margin,th1)
  
  fit <- mev::fit.gpd(
    margin,
    threshold = u,
    method = "Grimshaw",
    show = FALSE,
    MCMC = NULL,
    k = 4,
    tol = 1e-08,
    fpar = NULL,
    warnSE = FALSE
  )
  q0 <- evd::qgpd(1-1/(300*period*(1-th1)), 
                  scale = fit$estimate[1], 
                  shape =  fit$estimate[2] ) + u[1]
  q0 <- as.double(q0)
  
  fit <- fevd(margin, 
              threshold = u, 
              type = "GP", 
              verbose = TRUE)
  
  ## Not run: 
  ci <- ci(fit, type = "return.level", #method = "proflik",
           return.period = (300*period/365),
           xrange = estimates, verbose = plot)
  
  return(ci)
}

