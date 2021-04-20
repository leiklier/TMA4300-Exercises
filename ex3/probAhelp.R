## ARp.beta.est
##
## Estimate the parameters of an observed AR(p)-process, using least
## sum of squared residuals (LS) and least sum of absolute residuals
## (LA)
##
## The model is
##
##   x[t] = x[t-1]*beta[1] + ... + x[t-p]*beta[p] + e[t],  t=1,2,...,T
##
ARp.beta.est = function(x, p) {
    T = length(x)

    y = x[(1+p):T]
    C = matrix(x[rep(0:(T-p-1), p) + rep(p:1, each=T-p)], T-p, p)

    ## Least squares:
    beta.hat.LS = qr.solve(C,y)

    ## Least absolute deviations:
    if (p==1) {
        theta.hat = (optimize(function(theta,y,C) sum(abs(y-C%*%theta)),
                              c(-2,2), ## Don't require stationarity.
                              y=y, C=C))
        theta.hat = theta.hat$minimum
    } else {
        theta.hat = (optim(beta.hat.LS, ## Useful starting point.
                           function(theta,y,C) sum(abs(y-C%*%theta)),
                           gr=NULL,
                           y=y, C=C))
        theta.hat = theta.hat$par
    }
    beta.hat.LA = theta.hat

    return(list(LS=beta.hat.LS,
                LA=beta.hat.LA))
}

## ARp.filter
##
## Starting from the sequence in x0, calculate the AR(p)-sequence
##
##   x[t] = x[t-1]*beta[1] + ... + x[t-p]*beta[p] + e[t],  t=1,2,...,T
##
## where p = length(beta) and T = length(e)
##
## Hint: Use sample(e.observed, size=T, replace=TRUE) to generate a
##   resampled sequence of residuals to be used as input to this
##   function.
##
ARp.filter = function(x0, beta, e) {

    rx0 = rev(x0)
    if (is.null(x0)) {
        x0 = rep(0,length(beta))
    }
    # Caution: the function does not add the intial sequence to the 
    # time series. Further it requires the reverse intial sequence
    # as argument.
    x = filter(e, beta, method="recursive", init=rx0)

    return(c(x0, x))
}

## ARp.resid
##
## Calculate the observed residuals of an AR(p)-sequence,
##
##   e[t] = x[t] - x[t-1]*beta[1] - ... - x[t-p]*beta[p],  t=p+1,p+2,...,T
##
## where p = length(beta) and T = length(x)
##
## The output contains the values  e_{p+1} through e_T, ans so has length==T-p
##
ARp.resid = function(x, beta) {
    e.hat = filter(x,c(1,-beta),sides=1)[-(1:length(beta))]
    # re-center the residuals around zero
    return(e.hat - mean(e.hat))
}
