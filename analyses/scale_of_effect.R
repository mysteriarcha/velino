# Gaussian kernel analysis to find the "scale of effect" of altitude on
# vegetation diversity

# Our problem is to find which is the scale at which the effect of the
# predictor (altitude) on the response variable (alpha diversity) is
# strongest. We can address this problem with a (geographic) distance matrix 
# and a metric of evaluation of regression models.

# What model do we have to make? Let's check the 
lapply(sr_split, function(x) performance::check_distribution(x$alpha))
  
# Let's define a 
nll <- function(par, cov, y) {
  alpha <- par[1]
  beta <- par[2]
  lp <- alpha + beta*cov #linear predictor
  p <- plogis(lp) #back-transform
  loglike <- -sum(y*log(p) + (1-y)*log(1-p)) #negative jj
  return(loglike)
}
speciesrich$size <- as.numeric(speciesrich$size)
speciesrich$quota

lr.buffer <- optim(par = c(0, 0), fn = nll, 
                   cov = speciesrich$size,
                   y   = speciesrich$alpha, 
                   hessian = T
                   )

lr.buffer.vc <- solve(lr.buffer$hessian) #var-cov matrix
lr.buffer.se <- sqrt(diag(lr.buffer.vc)) #SE
lr.buffer.se



