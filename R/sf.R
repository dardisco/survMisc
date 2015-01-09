##' @name sf
##' @export sf
##'
##' @title Estimates of \bold{s}urvival (or hazard) \bold{f}unction
##' based on \eqn{n} and \eqn{e}
##'
##' @include tne.R
##' 
##' @param n \bold{N}umber at risk per time point (a \code{vector}).
##' @param e Number of \bold{e}vents per time point (a \code{vector}).
##' @param what What to return. See \bold{Value}.
##' 
##' @return
##' The return value will be a \code{vector},
##' unless \code{what="all"} (the default),
##' in which case it will be a \code{data.table}.
##' \cr
##' If \code{what="s"}, the \bold{s}urvival is returned, based on the
##' Kaplan-Meier or product-limit estimator.
##' This is \eqn{1} at \eqn{t=0} and thereafter is given by:
##' \deqn{\hat{S}(t) = \prod_{t \leq t_i} (1-\frac{e_i}{n_i} )}{
##'       S = PROD 1-e/n}
##'
##' If \code{what="sv"}, the \bold{s}urvival \bold{v}ariance is returned.
##' \cr
##' Greenwoods estimator of the variance of the Kaplan-Meier (product-limit)
##' estimator is:
##' \deqn{Var[\hat{S}(t)] = [\hat{S}(t)]^2 \sum_{t_i \leq t} \frac{e_i}{n_i (n_i-e_i)} }{
##'       Var[S](t) = S(t)^2 SUM (e) / n(n-e) }
##' 
##' If \code{what="h"}, the \bold{h}azard is returned, based on the the Nelson-Aalen estimator.
##' This has a value of \eqn{\hat{H}=0}{H=0} at \eqn{t=0} and thereafter is given by:
##' \deqn{\hat{H}(t) = \sum_{t \leq t_i} \frac{e_i}{n_i}  }{
##'       H(t) = SUM e/n}
##'
##' If \code{what="hv"}, the \bold{h}azard \bold{v}ariance is returned.
##' \cr
##' The variance of the Nelson-Aalen estimator is given by:
##' \deqn{Var[\hat{H}(t)] = \sum_{t_i \leq t} \frac{e_i}{n_i^2} }{
##'       SUM e/(n^2)}
##'
##' If \code{what="all"} (the default), \emph{all} of the above are returned
##' in a \code{data.table}, along with:
##' \cr
##' Survival, based on the Nelson-Aalen estimator. Given by
##'              \deqn{\hat{S_{na}} = \exp(-\hat{H})}{
##'                    S = exp(-H)}
##'           where \eqn{\hat{H}}{H} is hazard.
##' \cr \cr
##' HKM Hazard, based on the Kaplan-Meier estimator. Given by
##'              \deqn{\hat{H_{km}} = -\log \hat{S} }{
##'                    H = -log(S)}
##' where \eqn{\hat{S}}{S} is survival.
##'  
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' t1 <- tne(Surv(time=b1$t2, event=b1$d3))
##' sf(n=t1$n, e=t1$e, what="all")
##' ###
##' data(drug6mp, package="KMsurv")
##' s1 <- Surv(time=drug6mp$t2, event=drug6mp$relapse) # 6MP patients
##' t1 <- tne(s1)
##' sf(n=t1$n, e=t1$e, what="sv")
##' sf(n=t1$n, e=t1$e, what="hv")
##' 
##' @references Examples for are from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' \cr
##' \code{what="sv"}: Table 4.1A, pg 93.
##' \cr
##' \code{what="hv"}: Table 4.2, pg 94.
##' \cr
##' \code{what="all"}: Table 4.3, pg 97.
##' 
##'
sf <- function(n, e, what=c("all", "s", "sv", "h", "hv")){
    res1 <- s <- sv <- h <- hv <- sna <- hkm <- NULL
    if (!is.numeric(e)|!is.numeric(n)) stop("Vectors must be numeric")
    if (!length(e)==length(n)) stop("Vectors must be of equal length")
    what <- match.arg(what)
    if(what == "all"){
        res1 <- data.table::data.table(s = km(n, e))
        res1[, sv := kmv(n, e)]
        res1[, h := na(n, e)]
        res1[, hv := nav(n, e)]
        res1[, sna := exp(-h)]
        res1[, hkm := -log(s)]
        return(res1)
    }
    res1 <- switch(what,
                   s=km(n, e),
                   sv=kmv(n, e),
                   h=na(n, e),
                   hv=nav(n, e))
    return(res1)
}
###
### functions called above
### 
km <- function(n, e) cumprod(1 - (e / n))
kmv <- function(n, e) cumprod(1 - (e / n))^2 * cumsum( e / (n * (n - e)) ) 
na <- function(n, e) cumsum(e / n)
nav <- function(n, e) cumsum(e / (n^2))
