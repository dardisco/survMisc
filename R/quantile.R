##' @name quantile
##' @rdname quantile
##' @title Quantiles and median for \code{Surv}, \code{survfit} and \code{coxph} objects
##' @description Extends \code{stats::quantile} and \code{stats::quantile} to
##' work with \code{Surv}, \code{survfit} and \code{coxph} objects.
##' @export
##' 
quantile <- function(x, ...){
    UseMethod("quantile")
}
##'
##' @include sf.R
##' @include tne.R
##' 
##' @param x A \code{Surv}, \code{survfit} or \code{coxph} object.
##' @param ... Additional arguments (not implemented).
##' @param q (for \code{quantile}) Vector of quantiles
##' (expressed as percentage). For the \code{median}, \code{q=50}.
##' @param CI Include confidence interval.
##' \cr
##' Defaults are \code{CI=TRUE} for \code{quantile} and
##' \code{CI=FALSE} for \code{median}.
##' @param alpha Significance level \eqn{\alpha}{alpha}.
##' @param ci \bold{C}onfidence \bold{i}nterval.
##' \cr
##' One of: \bold{log} (the default), \bold{lin}ear or \bold{a}rcsine-\bold{s}quare \bold{r}oot.
##' @return For \code{quantile}:
##' A \code{data.table} (or a \code{list} of \code{data.table}s, one per stratum),
##' with columns:
##'   \item{q}{quantile}
##'   \item{t}{time}
##' If \code{CI = TRUE} then upper and lower confidence
##' intervals, as per argument \code{ci}).
##'  \item{l}{lower confidence limit}
##'  \item{u}{upper confidence limit}
##' For \code{median}:
##' A \code{data.table} with columns:
##'   \item{t}{time}
##'   \item{s}{stratum}
##' If \code{CI = TRUE} then a \code{list} of
##' \code{data.table}s, one per stratum, as above.
##'
##' @note If a time cannot be calculated, \code{NaN} is returned.
##' 
##' @seealso
##'
##' Confidence intervals are calculated as shown in the pointwise confidence intervals
##' in \code{\link{ci}}.
##'
##' @references Examples for quantiles are from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 4.2, pg 121.
##' 
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
##' 
##' @rdname quantile
##' @aliases quantile.Surv
##' @method quantile Surv
##' @export
quantile.Surv <- function(x, ...,
                          q=c(25, 50, 75),
                          CI=TRUE,
                          alpha=0.05,
                          ci=c("log", "lin", "asr")
                          ){
    if(!inherits(x, "Surv")) stop("Only applies to class 'Surv'")
    if(!attr(x,which="type")=="right") warning("Only applies to right censored data")
###
    dt1 <- tne(x)
    dt1[, s := sf(n=dt1[, n], e=dt1[, e], what="s")]
    dt1[, sv := sf(n=dt1[, n], e=dt1[, e], what="sv")]
### express as percentage
    p1 <- q / 100
### find time corresponding to quantile
    t1 <- sapply(p1, .findT, dt=dt1)
### convert alpha to z value
    z1 <- stats::qnorm(1 - alpha / 2)
### confidence interval
    ci <- match.arg(ci)
### get ranges for methods
    matZ1 <- as.matrix(switch(ci,
                              log = sapply(p1, .findLog, dt=dt1),
                              lin = sapply(p1, .findLin, dt=dt1),
                              asr = sapply(p1, .findArc, dt=dt1),
                              ), drop=FALSE)
    res1 <- data.table(q = q,
                       t = t1)
    if(CI){
        res1[, l := .findMin(matZ=matZ1, z=z1, dt=dt1)]
        res1[, u := .findMax(matZ=matZ1, z=z1, dt=dt1)]
    }
    setattr(res1, "ci", ci)
    return(res1)
}
##'
##' @rdname quantile
##' @aliases quantile.survfit
##' @method quantile survfit
##' @export
##' 
##' @examples
##' s1 <- survfit(Surv(t2, d3) ~ group, data=bmt)
##' quantile(s1)
##' 
quantile.survfit <-  function(x,
                              ...,
                              q=c(25,50,75),
                              CI=TRUE,
                              alpha=0.05,
                              ci=c("log", "lin", "asr")){
    if(!inherits(x, "survfit")) stop("Only applies to class 'survfit'")
    dt1 <- data.table(strata=summary(x)$strata,
                     t=summary(x)$time,
                     n=summary(x)$n.risk,
                     e=summary(x)$n.event,
                     s=summary(x)$surv)
### express as percentage
    p1 <- q / 100
### find time corresponding to quantile
    t1 <- sapply(p1, .findT, dt=dt1)
### convert alpha to z value
    z1 <- stats::qnorm(1 - alpha / 2)
### Greenwoods estimate of variance (from survEst)
    dt1[, sv := sf(n=dt1[, n], e=dt1[, e], what="sv")]
    ci <- match.arg(ci)
### strata
    s1 <- levels(dt1[, strata])
    res1 <- vector(mode="list", length=length(s1))
    names(res1) <- s1
    for (i in seq(s1)){
        res1[[i]] <- data.table(q = q,
                                t = t1)
        if(CI){
### ### ### get ranges for confidence interval
            matZ1 <- as.matrix(switch(ci,
                                      log = sapply(p1, .findLog, dt=dt1[strata==s1[i], ]),
                                      lin = sapply(p1, .findLin, dt=dt1[strata==s1[i], ]),
                                      asr = sapply(p1, .findArc, dt=dt1[strata==s1[i], ])
                                      ), drop=FALSE)
            res1[[i]][, l := .findMin(matZ=matZ1, z=z1,
                                  dt=dt1[strata==s1[i], ])]
            res1[[i]][, u := .findMax(matZ=matZ1, z=z1,
                            dt=dt1[strata==s1[i], ])]
        }
    }
    setattr(res1, "ci", ci)
    return(res1)
}
##'
##' @rdname quantile
##' @aliases quantile.coxph
##' @method quantile coxph
##' @export
##' 
##' @examples
##' c1 <- coxph(Surv(t2, d3)~ group, data=bmt)
##' quantile(c1)
##' 
quantile.coxph <- function(x, ...,
                           q=c(25,50,75),
                           CI=TRUE,
                           alpha=0.05,
                           ci=c("log", "lin", "asr")){
    stopifnot(inherits(x, "coxph")) 
    f1 <- deparse(x$call)
    f1 <- sub("coxph", "survfit", f1)
    s1 <- eval(parse(text=f1))
    quantile(s1, q=q, CI=CI, alpha=alpha, ci=ci)
}
##'
##'
##' @rdname quantile
##' @export 
##'
median <- function(x, ...){
    UseMethod("median")
}
##'
##' @rdname quantile
##' @aliases median.Surv
##' @method median Surv
##' @export
##' 
##' @examples
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' median(s1)
##' median(s1, CI=TRUE)
##' 
median.Surv <- function(x, ...,
                        CI=FALSE,
                        alpha=0.05,
                        ci=c("log", "lin", "asr")) {
    if(!inherits(x, "Surv")) stop("Only applies to Surv objects")
### for R CMD check
    s <- n <- e <- sv <- NULL
    ci <- match.arg(ci)
    quantile.Surv(x, q=50, alpha=alpha, ci=ci)
}
##'
##' @rdname quantile
##' @aliases median.survfit
##' @method median survfit
##' @export
##' 
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- survfit(Surv(t2, d3)~ group, data=bmt)
##' median(s1)
##' median(s1, ci="asr", CI=TRUE)
##'
median.survfit <- function(x, ...,
                           CI=FALSE,
                           alpha=0.05,
                           ci=c("log", "lin", "asr")) {
    if(!inherits(x, "survfit")) stop("Only applies to survfit objects")
    ci <- match.arg(ci)
    quantile.survfit(x, CI=CI, q=50, alpha=alpha, ci=ci)
    ## or can use shortcut here for most use cases
    ## read printout of survival:::print.survfit(x)
    ##     if (!CI) {
    ##         m1 <- read.table(textConnection(capture.output(x)),
    ##                          skip=2, header=TRUE)
    ##  get medians (drop -> prevent loss of rownames)
    ##         res1 <- data.table(t=m1[ ,"median"],
    ##                            s=rownames(m1))
    ##         return(res1)
    ##     }

}
##'
##' @rdname quantile
##' @aliases median.coxph
##' @method median coxph
##' @export
##' 
##' @examples
##' c1 <- coxph(Surv(t2, d3) ~ group, data=bmt)
##' median(c1)
median.coxph <- function(x, ...,
                         CI=FALSE,
                         alpha=0.05,
                         ci=c("log","lin","asr")) {
    f1 <- deparse(x$call)
    f1 <- sub("coxph", "survfit", f1)
    s1 <- eval(parse(text=f1))
    ci <- match.arg(ci)
    median.survfit(s1, CI=CI, alpha=alpha, ci=ci)
}
###
###----------------------------------------
###----------------------------------------
### functions used by quantile AND median
###----------------------------------------
###----------------------------------------
###
### in functions below:
### dt = data table with columns:
###   s = survival
###   sv = survival variance
###   t = time
### p = propbability
### s = survival
### 
### find time corresponding to quantile
### 
.findT <- function(p, dt) {
    ifelse(length(w1 <- which(dt[, s] <= p)),
           min(dt[w1, t]),
           NaN)
}
### 
### confidence intervals
### 
### linear interval
    .findLin <- function(p, dt) (dt[, s] - p) / sqrt((dt[, sv]))
### log transform (negative log log)
    .findLog <- function(dt, p){
        ( (log(-log(dt[, s])) - (log(-log(p))) ) * dt[, s] * log(dt[, s])) / sqrt(dt[, sv])
        }
### arcsin-square root transform
    .findArc <- function(p, dt){
        ( 2 * (asin(sqrt(dt[, s])) - asin(sqrt(p))) * sqrt( dt[, s] * (1-dt[, s])) ) / sqrt(dt[, sv])
    }
###
### in functions below:
### matZ = matrix of Z values, one column per quantile
### z = Z values to test (vector)
### dt = data table with columns
###   t = time
### 
### Find minimum time corresponding to z-value
### if any values < z
### get longest (largest) corresponding time < z
.findMin <- function(matZ, z, dt){
    sapply(seq(dim(matZ)[2]),
           function (i){
               ifelse(length(w1 <- min(which(matZ[, i] <= z))),
                      dt[w1, t],
                      NaN)
           })
}
### Find maximum time, as above
.findMax <- function(matZ, z, dt){
    sapply(seq(dim(matZ)[2]),
           function (i){
               ifelse(length(w1 <- which(matZ[, i] <= -z)),
                      dt[w1, t],
                      NaN)
           })
}
###----------------------------------------
###
### for R CMD check
    s <- n <- e <- sv <- u <- l <- NULL
