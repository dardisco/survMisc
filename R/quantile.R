##' @name quantile
##' @rdname quantile
##' @export quantile
quantile <- function(x, ...){
    UseMethod("quantile")
}
##' @rdname quantile
##' @include calcSurv.R
##' @include tne.R
##' @include survEst.R
##' @title Quantiles and median for \code{Surv} and \code{survfit} objects
##' @param x A \code{Surv} or \code{survfit} object
##' @param ... Additional arguments
##' @param q (for \code{quantile}) Vector of quantiles
##' (expressed as percentage).
##' @param alpha Significance level \eqn{\alpha}{alpha}
##' @param CI (for \code{median}) include confidence interval
##' @param method (for \code{median}) method for confidence interval, as below.
##' \cr
##' One of \code{log}, \code{linear} or \code{arcsine-square root}.
##' @return For \code{quantile}:
##' \cr
##' For a \code{Surv} object, a \code{data.frame} with quantile,
##' and upper and lower confidence intervals using the 3 methods.
##' For a \code{survfit} object, a list with
##' one element for each stratum. Each element is a
##' \code{data.frame} as above.
##' \cr \cr
##' For \code{median}:
##' \cr
##'  For a \code{Surv} object, the median time. For a \code{survfit} object
##' a \code{data.frame} with one row for each stratum.
##' \cr
##' If \code{CI} = \code{TRUE} then upper and lower confidence
##' intervals with one method (default is \code{log}).
##' @details
##'
##' Confidence intervals are calculated from \eqn{\sigma}{sigma} which is:
##' \deqn{\sigma (t) = \sqrt{ \frac{Var[\hat{S}(t)]}{\hat{S}^2(t)}} }{
##' ( Var[S(t)] / (S(t)^2) )^0.5}
##' The intervals are:
##' \cr
##'
##' \itemize{
##'
##' \item linear
##' \deqn{ \hat{S}(t) \pm Z_{1- \alpha} \sigma (t) \hat{S}(t)}{
##' S(t)+- Z(1-alpha) sigma(t) S(t)}
##' Where \eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival estimate.
##' \cr \cr
##'
##' \item log transform
##' \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
##' [S(t)^(1/theta), S(t)^theta]}
##' Where \eqn{\theta}{theta} is:
##' \deqn{ \exp{ \frac{Z_{1- \alpha} \sigma (t)}{ \log{\hat{S}(t)}}} }{
##' exp ( Z(1-alpha)sigma(t) / log(S(t)) )}
##' \cr \cr
##'
##' \item Arcsine-square root transform.
##' \cr
##' Upper:
##' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}} - \frac{Z_{1- \alpha}\sigma(t)}{2}
##' \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( max[0, arcsin S(t)^0.5 - Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##' Lower:
##' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}} +
##' \frac{Z_{1- \alpha}\sigma(t)}{2} \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##'  sin^2( min[pi/2, arcsin S(t)^0.5  + Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##'
##' }
##'
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' quantile(s1)
##' b1 <- bmt[bmt$group==2, ] # AML low-risk patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' quantile(s1)
##' b1 <- bmt[bmt$group==3, ] # AML high-risk patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' quantile(s1)
##' ###
##' @references Examples for quantiles are from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 4.2, pg 121.
##'
##' @export quantile.Surv
##' @aliases quantile.Surv
##' @method quantile Surv
##'
quantile.Surv <- function(x, ...,
                          q=c(25, 50, 75),
                          alpha=0.05){
    if(!inherits(x, "Surv")) stop("Only applies to class 'Surv'")
    if(!attr(x,which="type")=="right") warning("Applies to right censored data")
### calcSurv returns data.table by default
    s2 <- calcSurv(x)
    set(s2, j=6:9, value=NULL)
    class(s2) <- "data.frame"
    s2 <- s2[order(s2[ ,"t"]), ]
### express as percentage
    p1 <- q/100
### find time corresponding to quantile
    t1 <- sapply(p1, .findT, s=s2)
### convert alpha to z value
    z1 <- stats::qnorm(1-alpha/2)
### get ranges for methods
    Lin1 <- sapply(p1, .findLin, s=s2)
    Log1 <- sapply(p1, .findLog, s=s2)
    Arc1 <- sapply(p1, .findArc, s=s2)
### result
    res1 <- rbind(t1,
                .findMin(Lin1, z1, s=s2),
                .findMax(Lin1, z1, s=s2),
                .findMin(Log1, z1, s=s2),
                .findMax(Log1, z1, s=s2),
                .findMin(Arc1, z1, s=s2),
                .findMax(Arc1, z1, s=s2)
                )
    colnames(res1) <- q
    rownames(res1) <- c("quantile",
                      "lin.lower","lin.upper",
                      "log.lower","log.upper",
                      "a.s.lower","a.s.upper")
    return(res1)
}
##'
##' @rdname quantile
##' @export quantile.survfit
##' @aliases quantile.survfit
##' @method quantile survfit
##' @examples
##' s1 <- survfit(Surv(t2, d3)~ group, data=bmt)
##' quantile(s1)
##' ###
quantile.survfit <-  function(x, ..., q=c(25,50,75), alpha=0.05){
    if(!inherits(x, "survfit")) stop("Only applies to class 'survfit'")
    s1 <- data.frame(strata=summary(x)$strata,
                     t=summary(x)$time,
                     n=summary(x)$n.risk,
                     e=summary(x)$n.event,
                     SKM=summary(x)$surv)
### express as percentage
    p1 <- q/100
### find time corresponding to quantile
    t1 <- sapply(p1, .findT, s=s1)
### convert alpha to z value
    z1 <- stats::qnorm(1-alpha/2)
### Greenwoods estimate of variance (from survEst)
    s1$SKMV <- Gw(n=s1$n, e=s1$e)
### declare variables
    res <- vector(mode="list", length=length(levels(s1$strata)))
    names(res) <- levels(s1$strata)
    for (i in 1:length(levels(s1$strata))){
### subset per strata
        s2 <- s1[s1$strata==levels(s1$strata)[i], ]
### get ranges for methods
        Lin1 <- sapply(p1, .findLin, s=s2)
        Log1 <- sapply(p1, .findLog, s=s2)
        Arc1 <- sapply(p1, .findArc, s=s2)
        res1 <- rbind(t1,
                      .findMin(Lin1, z1, s=s2),
                      .findMax(Lin1, z1, s=s2),
                      .findMin(Log1, z1, s=s2),
                      .findMax(Log1, z1, s=s2),
                      .findMin(Arc1, z1, s=s2),
                      .findMax(Arc1, z1, s=s2)
                      )
        colnames(res1) <- q
        rownames(res1) <- c("quantile",
                            "lin.lower","lin.upper",
                            "log.lower","log.upper",
                            "a.s.lower","a.s.upper")
        res[[i]] <- res1
    }
    return(res)
}
##' @rdname quantile
##' @aliases quantile.coxph
##' @method quantile coxph
##' @S3method quantile coxph
##' @examples
##' c1 <- coxph(Surv(t2, d3)~ group, data=bmt)
##' quantile(c1)
quantile.coxph <- function(x, ..., q=c(25,50,75), alpha=0.05){
    f1 <- deparse(x$call)
    f1 <- sub("coxph", "survfit", f1)
    s1 <- eval(parse(text=f1))
    quantile(s1, q=q, alpha=alpha)
}
###
###----------------------------------------
### functions used by quantile AND median
###
### find time corresponding to quantile
    .findT <- function(p, s) {
        ifelse(any(s[, "SKM"] <= p),
               min(s[which(s[, "SKM"] <= p), "t"]),
               NA)
        }
### linear
    .findLin <- function(p, s) {
        (s[ ,"SKM"]-p) / sqrt((s[ ,"SKMV"]))
            }
### log transform (negative log log)
    .findLog <- function(p, s){
        ( ( log(-log(s[, "SKM"])) - (log(-log(p))) )*s[, "SKM"]*log(s[, "SKM"]))/ sqrt(s[, "SKMV"])
        }
### arcsin-square root transform
    .findArc <- function(p, s){
        ( 2*(asin(sqrt(s[, "SKM"])) - asin(sqrt(p)))*sqrt( s[, "SKM"]*(1-s[, "SKM"])) )/sqrt(s[, "SKMV"])
    }
### Find minimum time
    .findMin <- function(mat, z1, s){
        sapply(1:ncol(mat),
               function (i){
### if any values < z
                   ifelse(any(mat[ ,i]<=z1),
### get longest (largest) corresponding time < z
### from Surv object (in survfit = 1 strata)
                          s[min(which(mat[ ,i]<=z1)),"t"],
                          NA)
               }
               )
    }
### find maximum time
    .findMax <- function(mat, z1, s){
        sapply(1:ncol(mat),
               function (i){
                   ifelse(any(mat[ ,i]<=-z1),
                          s[max(which(mat[ ,i]<=-z1)),"t"],
                          NA)
               }
               )
    }
##'
##' @rdname quantile
##' @export median
##'
median <- function(x, ...){
    UseMethod("median")
}
##' @rdname quantile
##' @aliases median.Surv
##' @method median Surv
##' @S3method median Surv
##' @examples
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' median(s1)
##' median(s1, CI=TRUE)
median.Surv <- function(x, ...,
                        CI=FALSE,
                        alpha=0.05,
                        method=c("log", "lin", "asr")) {
    if(!inherits(x, "Surv")) stop("Only applies to Surv objects")
### calcSurv returns data.table by default
    s2 <- calcSurv(x)
    set(s2, j=6:9, value=NULL)
    class(s2) <- "data.frame"
    s2 <- s2[order(s2[, "t"]),]
### find time corresponding to 50th percentile
    t1 <- .findT(0.5, s2)
### get medians (drop -> prevent loss of rownames)
    if (!CI) return(t1)
###
### character match
    method <- match.arg(method)
### convert alpha to z value
    z1 <- stats::qnorm(1-alpha/2)
### CI based on method
    range1 <- as.matrix(switch(method,
                               log = .findLog(0.5, s=s2),
                               lin = .findLin(0.5, s=s2),
                               asr = .findArc(0.5, s=s2)
                               ), drop=FALSE)
    lower <- .findMin(range1, z1, s=s2)
    upper <- .findMax(range1, z1, s=s2)
    res1 <- cbind(t1, lower, upper)
    colnames(res1)[1] <- "median"
    return(res1)
}
##'
##' @rdname quantile
##' @aliases median.survfit
##' @method median survfit
##' @S3method median survfit
##'
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- survfit(Surv(t2, d3)~ group, data=bmt)
##' median(s1)
##' median(s1, CI=TRUE, method="asr")
##'
median.survfit <- function(x, ...,
                           CI=FALSE,
                           alpha=0.05,
                           method=c("log", "lin", "asr")) {
    if(!inherits(x, "survfit")) stop("Only applies to survfit objects")
### get printout of survival:::print.survfit(x)
    m1 <- read.table(textConnection(capture.output(x)),
                     skip=2, header=TRUE)
### get medians (drop -> prevent loss of rownames)
    m1 <- m1[ ,"median" , drop=FALSE]
    if (!CI) return(m1)
###
### character match
    method <- match.arg(method)
### convert alpha to z value
    z1 <- stats::qnorm(1-alpha/2)
### get estimates of survival at those times
    m1$strata <- rownames(m1)
    s1 <- data.frame(strata=summary(x)$strata,
                     t=summary(x)$time,
                     n=summary(x)$n.risk,
                     e=summary(x)$n.event,
                     SKM=summary(x)$surv)
### Greenwoods estimate of variance (from survEst)
    s1$SKMV <- Gw(n=s1$n, e=s1$e)
### declare variables
    lower <- upper <- NULL
    for (i in 1: length(m1$strata)){
### subset per strata
        s2 <- s1[s1$strata==m1$strata[i], ]
### CI based on method
        range1 <- as.matrix(switch(method,
                         lin = .findLin(0.5, s=s2),
                         log = .findLog(0.5, s=s2),
                         asr = .findArc(0.5, s=s2)
                         ), drop=FALSE)
        lower[i] <- .findMin(range1, z1, s=s2)
        upper[i] <- .findMax(range1, z1, s=s2)
    }
    res1 <- cbind(m1[, -2, drop=FALSE], lower, upper)
    return(res1)
}
##' @rdname quantile
##' @aliases median.coxph
##' @method median coxph
##' @S3method median coxph
##' @examples
##' c1 <- coxph(Surv(t2, d3) ~ group, data=bmt)
##' median(c1)
median.coxph <- function(x, ..., CI=FALSE,
                           alpha=0.05,
                           method=c("log","lin","asr")) {
    f1 <- deparse(x$call)
    f1 <- sub("coxph", "survfit", f1)
    s1 <- eval(parse(text=f1))
    median(s1, CI=CI, alpha=alpha, method=method)
}
