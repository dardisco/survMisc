##' @name ci
##' @title Confidence intervals for survival curves
##' @rdname ci
##' @export ci
##'
ci <- function(x, ...){
    UseMethod("ci")
    }
##' @rdname ci
##' @aliases ci.survfit
##' @method ci survfit
##' @S3method ci survfit
##' @include tne.R
##' @include calcSurv.R
##' @include survEst.R
##'
##' @param x An object of class \code{survfit}
##' @param ... Additional arguments (not implemented)
##' @param CI Confidence intervals. As the function currently relies on lookup
##' tables, currently only 95\% (the default), 90\% and 99\% are supported.
##' @param how Method to use for confidence interval. \code{point}
##' (the default) uses pointwise confirence intervals. The alternatives use
##' confidence \emph{bands} (see details).
##' @param trans Transformation to use. Default \code{lin} is linear. Also
##' supported are log and arcsine-square root transformations.
##' @param tL Lower time point. Used in construction of confidence bands.
##' @param tU Upper time point. Used in construction of confidence bands.
##' @return A \code{survfit} object. The \code{upper} and \code{lower}
##' elements in the list (representing confidence intervals)
##' are modified from the original. Other elements
##' will also be shortened if the time range under consideration has been
##' reduced from the original.
##'
##' @details
##' In the equations below
##' \deqn{\sigma^2_s(t) = \frac{\hat{V}[\hat{S}(t)]}{\hat{S}^2(t)} }{
##'  sigma^2(t) = V[S(t)]/[S(t)]^2}
##' Where \eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival estimate and
##' \eqn{\hat{V}[\hat{S}(t)]}{V[S(t)]} is Greenwood's estimate of its
##' variance.
##' \cr
##' The pointwise confidence intervals are valid for individual
##' times, e.g. median and \code{\link{quantile}} values.
##' When plotted and joined for multiple points they tend to
##' be narrower than the \emph{bands} described below.
##' Thus they tend to exaggerate the impression of certainty
##' when used to plot confidence intervals for a time range.
##' They should not be interpreted as giving the intervals
##' within which the \emph{entire} survival function lies.
##' \cr
##'
##' For a given significance level \eqn{\alpha}{alpha},
##' they are calculated using the standard normal distribution \eqn{Z}
##' as follows:
##'
##' \itemize{
##'
##' \item linear
##' \deqn{ \hat{S}(t) \pm Z_{1- \alpha} \sigma (t) \hat{S}(t)}{
##' S(t)+- Z(1-alpha) sigma(t) S(t)}
##'
##' \item log transform
##' \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
##' [S(t)^(1/theta), S(t)^theta]}
##' where
##' \deqn{ \theta = \exp{ \frac{Z_{1- \alpha} \sigma (t)}{ \log{\hat{S}(t)}}} }{
##' theta = exp ( Z(1-alpha)sigma(t) / log(S(t)) )}
##'
##' \item arcsine-square root transform.
##' \cr
##' Upper:
##' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}} -
##'   \frac{Z_{1- \alpha}\sigma(t)}{2}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2(max[0, arcsin S(t)^0.5 - Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##' Lower:
##' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}} +
##'   \frac{Z_{1- \alpha}\sigma(t)}{2}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2(min[pi/2, arcsin S(t)^0.5  + Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##'
##' }
##'
##' Confidence \emph{bands} give the values within which the survival function
##' falls within a \emph{range} of timepoints.
##' \cr \cr
##' The time range under consideration is given so that
##' \eqn{t_l \geq t_{min}}{tL >= min(t)}, the lowest event time and
##' \eqn{t_u \leq t_{max}}{tU <= max(t)}, the largest event time.
##' \cr
##' For a sample size \eqn{n} and \eqn{0 < a_l < a_u <1}:
##' \deqn{a_l = \frac{n\sigma^2_s(t_l)}{1+n\sigma^2_s(t_l)}}{
##' a_l = n*sigma^2(t_l) / [1+n*sigma^2(t_l)]}
##' \deqn{a_u = \frac{n\sigma^2_s(t_u)}{1+n\sigma^2_s(t_u)}}{
##' a_u = n*sigma^2(t_u) / [1+n*sigma^2(t_u)]}
##'
##' For the \bold{Nair} or \bold{equal precision} (\bold{EP}) confidence bands,
##' we begin by obtaining the relevant
##' confidence coefficient \eqn{c_{\alpha}}{c[alpha]}. This is obtained from
##' the upper \eqn{\alpha}{a}-th fractile of the random variable
##' \deqn{U = \sup{|W^o(x)\sqrt{[x(1-x)]}|, \quad a_l \leq x \leq a_u}}{
##' U = sup{ |W(x)[x(1-x)]^0.5|, a_l <= x <= a_u} }
##' Where \eqn{W^o}{W} is a standard Brownian bridge.
##' \cr
##' The intervals are:
##'
##' \itemize{
##'
##' \item linear
##' \deqn{ \hat{S}(t) \pm c_{\alpha} \sigma_s(t) \hat{S}(t)}{
##' S(t)+- c[alpha] sigma(t) S(t)}
##'
##' \item log transform
##' \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
##' [S(t)^(1/theta), S(t)^theta]}
##' where
##' \deqn{ \theta = \exp{ \frac{c_{\alpha} \sigma_s(t)}{ \log{\hat{S}(t)}}} }{
##' theta = exp ( c[alpha]*sigma(t) / log(S(t)) )}
##'
##' \item arcsine-square root transform.
##' \cr
##' Upper:
##' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}}
##'   - \frac{c_{\alpha}\sigma_s(t)}{2}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( max[0, arcsin S(t)^0.5 - c[alpha]*sigma(t)/2 (S(t)/1-S(t))^0.5])}
##' Lower:
##' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}}
##'   + \frac{c_{\alpha}\sigma_s(t)}{2}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( min[pi/2, arcsin S(t)^0.5 - c[alpha]*sigma(t)/2 (S(t)/1-S(t))^0.5])}
##'
##' }
##'
##' For the \bold{Hall-Wellner} bands the confidence coefficient
##' \eqn{k_{\alpha}}{k[alpha]}
##' is obtained from the upper \eqn{\alpha}{a}-th fractile of a
##' Brownian bridge.
##' \cr
##' In this case \eqn{t_l} can be \eqn{=0}.
##' The intervals are:
##'
##' \itemize{
##'
##' \item linear
##' \deqn{ \hat{S}(t) \pm
##' k_{\alpha} \frac{1+n\sigma^2_s(t)}{\sqrt{n}} \hat{S}(t)}{
##' S(t)+- k[alpha] [1+n*sigma^2(t)]*S(t) / n^0.5 }
##'
##' \item log transform
##' \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
##' [S(t)^(1/theta), S(t)^theta]}
##' where
##' \deqn{ \theta = \exp{ \frac{k_{\alpha}[1+n\sigma^2_s(t)]}{
##'    \sqrt{n}\log{\hat{S}(t)}}} }{
##' theta = exp ( k[alpha]*[1+n*sigma^2(t)] / n^0.5 * log(S(t)) )}
##'
##' \item arcsine-square root transform.
##' Upper:
##' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}}
##'   - \frac{k_{\alpha}[1+n\sigma_s(t)]}{2\sqrt{n}}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( max[0, arcsin S(t)^0.5 - k[alpha]*[1+n*sigma^2(t)]/(2*n^0.5) (S(t)/1-S(t))^0.5])}
##' Lower:
##' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}}
##'   + \frac{k_{\alpha}[1+n\sigma^2_s(t)]}{2\sqrt{n}}
##'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( min[pi/2, arcsin S(t)^0.5 - k[alpha]*[1+n*sigma^2(t)]/(2*n^0.5) (S(t)/1-S(t))^0.5])}
##'
##' }
##'
##' @note
##' The function is loosely based on \code{km.ci::km.ci}.
##' \cr
##' For the Nair and Hall-Wellner bands, the function currently
##' relies on the lookup tables in \code{package:km.ci}.
##' \cr
##' Generally, the arcsin-square root transform has the best coverage properties.
##' \cr
##' All bands have good coverage properties for samples as small as \eqn{n=20},
##' except for the \bold{Nair / EP} bands with a linear transformation,
##' which perform poorly when \eqn{n < 200}.
##'
##' @references
##' Nair V, 1984. Confidence bands for survival functions with censored data:
##' a comparative study. \emph{Technometrics}. \bold{26}(3):265-75.
##' \href{http://www.jstor.org/stable/1267553}{JSTOR}.
##' @references
##' Hall WJ, Wellner JA, 1980. Confidence bands for a survival curve
##' from censored data.
##' \emph{Biometrika}. \bold{67}(1):133-43.
##' \href{http://www.jstor.org/stable/2335326}{JSTOR}.
##' @references
##' Examples are from:
##' Klein J, Moeschberger M 2003.
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Section 4.4, pg 111.
##' @seealso
##' \code{\link{KM}}
##' \code{\link{quantile}}
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- survfit(Surv(t2, d3) ~ 1, data=bmt[bmt$group==1, ])
##' ci(s1, how="nair", trans="lin", tL=100, tU=600)
##' s2 <- survfit(Surv(t2, d3) ~ group, data=bmt)
##' ci(s2, CI="0.99", how="point", trans="asin", tL=100, tU=600)
ci.survfit <- function(x, ...,
                       CI=c("0.95", "0.9", "0.99"),
                       how=c("point", "nair", "hall"),
                       trans=c("lin", "log", "asin"),
                       tL=NULL,
                       tU=NULL){
    lower <- upper <- cA <- kA <- n <- e <- NULL
    SKM <- SKMV <- SKMse <- sigSq <- NULL
    stopifnot(inherits(x, "survfit"))
    trans <- match.arg(trans)
    CI <- 100*as.numeric(match.arg(CI))
    how <- match.arg(how)
###
    t1 <- tne(x, return="list", onlyEvents=FALSE)
    for(i in seq_along(t1)){
        t1[[i]] <- Surv(t1[[i]]$t, t1[[i]]$e)
    }
    t1 <- lapply(t1, calcSurv)
    n1 <- unlist(lapply(t1, function(x) max(x[, n])))
    lapply(t1, function(x) set(x, j=6L, value=x[, SKMV / SKM^2]))
    lapply(t1, function(x) set(x, j=7L:ncol(x), value=NULL))
    lapply(t1, function(x) setnames(x, c("t", "n", "e",
                                         "SKM", "SKMV", "sigSq")))
    if(is.null(tL)) tL <- min(unlist(lapply(t1, function(x) x[, t])))
    if(is.null(tU)) tU <- max(unlist(lapply(t1, function(x) x[, t])))
    t1 <- lapply(t1, function(x) x[t > tL & t < tU, j=1L:ncol(x), with=FALSE])
###
    linLow <- switch(how,
                     point = function(x, z) {
                         x[, SKM] - (z * sqrt(x[, sigSq]) * x[, SKM])
                     },
                     nair = function(x, cA) {
                         x[, SKM] - (cA * sqrt(x[, sigSq]) * x[, SKM])
                     },
                     hall = function(x, kA, n) {
                         x[, SKM] - (((kA * (1 + n * x[, sigSq]) / sqrt(n))) * x[, SKM])
                     })
    linUp <- switch(how,
                    point = function(x, z){
                        x[, SKM] + (z * sqrt(x[, sigSq]) * x[, SKM])
                    },
                    nair =  function(x, cA){
                        x[, SKM] + (cA * sqrt(x[, sigSq]) * x[, SKM])
                    },
                    hall =  function(x, kA, n){
                        x[, SKM] + (((kA * (1 + n * x[, sigSq]) / sqrt(n))) * x[, SKM])
                    })
    logLow <- function(x, theta){
        x[, SKM]^(1/theta)
    }
    logUp <- function(x, theta){
        x[, SKM]^(theta)
    }
    asinLow <- switch(how,
                      point =  function(x, z){
                          sin(asin(x[, sqrt(SKM)]) - (0.5 * z * x[, sqrt(sigSq)] *
                                                      sqrt(x[, SKM] / (1 - x[, SKM]))))^2
                      },
                      nair =  function(x, cA){
                          sin(asin(x[, sqrt(SKM)]) - (0.5 * cA * x[, sqrt(sigSq)] *
                                                      sqrt(x[, SKM] / (1 - x[, SKM]))))^2
                      },
                      hall =  function(x, kA, n){
                          sin(asin(x[, sqrt(SKM)]) -
                              (0.5 * kA * (1 + n * x[, sigSq]) / sqrt(n)) *
                              sqrt(x[, SKM] / (1 - x[, SKM])))^2
                      })
    asinUp <- switch(how,
                     point =  function(x, z){
                         sin(asin(x[, sqrt(SKM)]) + (0.5 * z * x[, sqrt(sigSq)] *
                                                     sqrt(x[, SKM] / (1 - x[, SKM]))))^2
                     },
                     nair = function(x, cA){
                         sin(asin(x[, sqrt(SKM)]) + (0.5 * cA * x[, sqrt(sigSq)] *
                                                     sqrt(x[, SKM] / (1 - x[, SKM]))))^2
                     },
                     hall =  function(x, kA, n){
                         sin(asin(x[, sqrt(SKM)]) +
                             (0.5 * kA * (1 + n * x[, sigSq]) / sqrt(n)) *
                             sqrt(x[, SKM] / (1 - x[, SKM])))^2
                     })
     if(how=="point"){
### get Z (quantile from normal distribution)
         alpha <- (100 - CI) / 100
         z1 <- stats::qnorm(1 - alpha/2)
         if(trans=="log"){
            theta1 <- vector(mode="list", length=length(n1))
            for (i in seq_along(n1)){
                theta1[[i]] <- exp(z1 * sqrt(t1[[i]][, sigSq]) /
                                   log(t1[[i]][, SKM]))
            }
        }
        for(i in seq_along(n1)){
            t1[[i]][, lower :=
                    switch(trans,
                           lin = linLow(x=t1[[i]], z=z1[1]),
                           log = logLow(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinLow(x=t1[[i]], z=z1[1])
                           )]
            t1[[i]][, lower := ifelse(lower < 0, 0, lower)]
            t1[[i]][, upper :=
                    switch(trans,
                           lin = linUp (x=t1[[i]], z=z1[1]),
                           log = logUp(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinUp(x=t1[[i]], z=z1[1])
                           )]
            t1[[i]][, upper := ifelse(upper > 1, 1, upper)]
        }
     }
    if(how=="nair" | how=="hall"){
### get lookup table for confidence coefficient
        d1s <- paste0("critical.value.", how, ".", CI)
        do.call(data, list(eval(substitute(d1s)), package="km.ci"))
        do.call(assign, list("d1", eval(parse(text=d1s))))
### generate a, lower and upper
        genA <- function(n, sigSq) round((n * sigSq) / (1 + n * sigSq), 1)
        aL <- aU <- vector(mode="double", length=length(n1))
        for (i in seq_along(n1)){
            aL[i] <- genA(n=n1[i], sigSq=t1[[i]][1, sigSq])
            aU[i] <- genA(n=n1[i], sigSq=t1[[i]][nrow(t1[[i]]), sigSq])
        }
    }
###
    if(how=="nair"){
        rownames(d1) <- seq(0.1, 0.98, by=0.02)
        colnames(d1) <- seq(0.02, 0.6, by=0.02)
###
        err1 <- "a_L &/ a_U are outside the range of lookup table. \n Suggest try narrower range for time i.e. increase tL &/ decrease tU"
        if(!all(c(aU %in% rownames(d1), aL %in% colnames(d1)))) stop(err1)
### get confidence coefficient
        cA1 <- diag(d1[rownames=as.character(aU),
                       colnames=as.character(aL)])
        if(trans=="log"){
            theta1 <- vector(mode="list", length=length(n1))
            for (i in seq_along(n1)){
                theta1[[i]] <- exp(cA[i] * sqrt(t1[[i]][, sigSq]) /
                                   log(t1[[i]][, SKM]))
            }
        }
        for(i in seq_along(n1)){
            t1[[i]][, lower :=
                    switch(trans,
                           lin = linLow(x=t1[[i]], cA=cA1[i]),
                           log = logLow(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinLow(x=t1[[i]], cA=cA1[i])
                           )]
            t1[[i]][, lower := ifelse(lower < 0, 0, lower)]
            t1[[i]][, upper :=
                    switch(trans,
                           lin = linUp(x=t1[[i]], cA=cA1[i]),
                           log = logUp(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinUp(x=t1[[i]], cA=cA1[i])
                           )]
            t1[[i]][, upper := ifelse(upper > 1, 1, upper)]
        }
    }
###
    if(how=="hall"){
        rownames(d1) <- seq(0.1, 1.0, by=0.02)
        colnames(d1) <- seq(0, 0.6, by=0.02)
        err1 <- "a_L &/ a_U are outside the range of lookup table. \n Suggest try narrower range for time i.e. increase tL &/ decrease tU"
        if(!all(c(aU %in% rownames(d1), aL %in% colnames(d1)))) stop(err1)
### get confidence coefficient
        kA1 <- diag(d1[rownames=as.character(aU), colnames=as.character(aL)])
        if(trans=="log"){
            theta1 <- vector(mode="list", length=length(n1))
            for (i in seq_along(n1)){
                theta1[[i]] <- exp(kA[i] * (1 + n1[i] * t1[[i]][, sigSq])
                                   / (sqrt(n1[i]) * log(t1[[i]][, SKM])))
            }
        }
        for(i in seq_along(n1)){
            t1[[i]][, lower :=
                    switch(trans,
                           lin = linLow(x=t1[[i]], kA=kA1[i], n=n1[i]),
                           log = logLow(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinLow(x=t1[[i]], kA=kA1[i], n=n1[i])
                           )]
            t1[[i]][, lower := ifelse(lower < 0, 0, lower)]
            t1[[i]][, upper :=
                    switch(trans,
                           lin = linUp(x=t1[[i]], kA=kA1[i], n=n1[i]),
                           log = logUp(x=t1[[i]], theta=theta1[[i]]),
                           asin = asinUp(x=t1[[i]], kA=kA1[i], n=n1[i])                           )]
            t1[[i]][, upper := ifelse(upper > 1, 1, upper)]
        }
    }
    for(i in seq_along(t1)){
### add column for censored
        t1[[i]][, "c" := c( - diff(n) - e[1:(length(e) - 1)],
                      !sign(e)[length(e)])]
    }
### repackage as survfit object
    res <- vector("list", 14)
    names(res) <- c("n", "time", "n.risk", "n.event", "n.censor",
                    "surv",
                    "type", "strata",
                    "std.err", "upper", "lower",
                    "conf.type", "conf.int", "call")
    res[[1]] <- n1
    res[[2]] <- unname(unlist(lapply(t1, function(x) x[, t])))
    res[[3]] <- unname(unlist(lapply(t1, function(x) x[, n])))
    res[[4]] <- unname(unlist(lapply(t1, function(x) x[, e])))
    res[[5]] <- unname(unlist(lapply(t1, function(x) x[, c])))
    res[[6]] <- unname(unlist(lapply(t1, function(x) x[, SKM])))
    res[[7]] <- "right"
    res[[8]] <- unlist(lapply(t1, function(x) nrow(x)))
### standard error is square root of variance
    res[[9]] <- unname(unlist(lapply(t1, function(x) x[, sqrt(SKMV)])))
    res[[10]] <- unname(unlist(lapply(t1, function(x) x[, upper])))
    res[[11]] <- unname(unlist(lapply(t1, function(x) x[, lower])))
    res[[12]] <- paste(trans, how, sep="_")
    res[[13]] <- CI/100
    res[[14]] <- x$call
    class(res) <- "survfit"
    return(res)
}
