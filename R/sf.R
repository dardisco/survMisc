#' @name sf
#' @export sf
#'
#' @title Estimates of \bold{s}urvival (or hazard) \bold{f}unction
#' based on \eqn{n} and \eqn{e}
#'
#' @include tne.R
#'
#' @param x
#' @param ... Not implemented
#' @param n \bold{N}umber at risk per time point (a \code{vector})
#' @param what See return, below
#'
#' @return
#' The return value will be a \code{vector},
#' unless \code{what="all"} (the default),
#' in which case it will be a \code{data.table}.
#' \cr
#' If \code{what="s"}, the \bold{s}urvival is returned, based on the
#' Kaplan-Meier or product-limit estimator.
#' This is \eqn{1} at \eqn{t=0} and thereafter is given by:
#' \deqn{\hat{S}(t) = \prod_{t \leq t_i} (1-\frac{e_i}{n_i} )}{
#'       S[t] = prod (1 - e[t]) / n[t] }
#'
#' If \code{what="sv"}, the \bold{s}urvival \bold{v}ariance is returned.
#' \cr
#' Greenwoods estimtor of the variance of the Kaplan-Meier (product-limit)
#' estimator is:
#' \deqn{Var[\hat{S}(t)] = [\hat{S}(t)]^2 \sum_{t_i \leq t} \frac{e_i}{n_i (n_i-e_i)} }{
#'       Var(S)[t] = S[t]^2 sum e[t] / (n[t] * (n[t] - e[t])) }
#'
#' If \code{what="h"}, the \bold{h}azard is returned, based on the the Nelson-Aalen estimator.
#' This has a value of \eqn{\hat{H}=0}{H=0} at \eqn{t=0} and thereafter is given by:
#' \deqn{\hat{H}(t) = \sum_{t \leq t_i} \frac{e_i}{n_i}  }{
#'       H[t] = sum e[t] / n[t]}
#'
#' If \code{what="hv"}, the \bold{h}azard \bold{v}ariance is returned.
#' \cr
#' The variance of the Nelson-Aalen estimator is given by:
#' \deqn{Var[\hat{H}(t)] = \sum_{t_i \leq t} \frac{e_i}{n_i^2} }{
#'       SUM e/(n^2)}
#'
#' If \code{what="all"} (the default), \emph{all} of the above are returned
#' in a \code{data.table}, along with:
#' \cr
#' Survival, based on the Nelson-Aalen estimator. Given by:
#'  \deqn{\hat{S_{na}}=e^{H}}{
#'        S[t] = exp H[t]}
#' where \eqn{H} is hazard.
#' HKM Hazard, based on the Kaplan-Meier estimator. Given by:
#'  \deqn{\hat{H_{km}}=-\log{S}}{
#'        H[t] = -log S[t]}
#' where \eqn{S} is survival.
#'
#' @examples
#' ## K&M 2nd ed. Table 4.1A, pg 93.
#' ## 6MP patients
#' data("drug6mp", package="KMsurv")
#' s1 <- with(drug6mp, Surv(time=t2, event=relapse))
#' t1 <- tne(s1)
#' sf(x=t1$e, n=t1$n, what="sv")
#' ## K&M 2nd ed. Table 4.2, pg 94.
#' data("bmt", package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' t2 <- tne(Surv(time=b1$t2, event=b1$d3))
#' with(t2, sf(x=e, n=n, what="Hv"))
#' ## K&M 2nd ed. Table 4.3, pg 97.
#' sf(x=t2$e, n=t2$n, what="all")
#'
sf <- function(x, ...) UseMethod("sf")
###
sf.default <- function(x, ...){
    stopifnot(all(x >= 0 && x <=1))
    t1 <- tne(x)
    return(sf.tne(t1))
}
#' @examples
#' t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney)
#' sf(t1)
#' t2 <- tne(t1)
#' sf(t2)[order(cg), ]
sf.tne <- function(x, ...,
                   what=c("S", "H"),
                   sigSq=FALSE){
    stopifnot(inherits(x, "tne"))
    what <- match.arg(what)
    ## functions to use
    if (what=="S"){
        fun1 <- km
        fun1v <- kmv
    } else {
        fun1 <- na
        fun1v <- nav
    }
    ## name of variance
    nv1 <- paste0(what, "v")
    if (attr(x, "shape")=="long"){
        if (attr(x, "ncg") < 2){
            res1 <- data.table::data.table(
                x[, t],
                x[, fun1(e, n)],
                x[, fun1v(e, n)])
            data.table::setnames(res1, c("t", what, nv1))
            if (sigSq) res1[, "sigSq" := Sv / S^2]
        } else {
            res1 <- data.table::data.table(
                x[, t],
                x[, cg],
                x[, fun1(e, ncg), by=cg][, V1],
                x[, fun1v(e, ncg), by=cg][, V1])
        data.table::setnames(res1, c("t", "cg", what, nv1))
        if (sigSq) res1[, "sigSq" := Sv / S^2, by=cg]
        }
    }
    if(attr(x, "shape")=="wide"){
        e_ <- grep("e_", names(x))
        n_ <- grep("n_", names(x))
        res1 <- data.table::data.table(
            x[, t],
            mapply(FUN=fun1, x[, .SD, .SDcols=e_], x[, .SD, .SDcols=n_]),
            mapply(FUN=fun1v, x[, .SD, .SDcols=e_], x[, .SD, .SDcols=n_]))
        ## names of covariate groups
        if (attr(x, "abbNames")){
            n1 <- attr(x, "longNames")[, id]
        } else {
            n1 <- attr(x, "longNames")[, longName]
        }
        data.table::setnames(res1,
                             c("t",
                               as.vector(t(outer(c(what, nv1), n1, paste0)))))
        if (sigSq){
            ssN1 <- paste0("sigSq", n1) 
            res1[, (ssN1) := mapply(function(x, y) list(x / y^2),
                              res1[, .SD, .SDcols=grep("Sv", names(res1))],
                              res1[, .SD, .SDcols=grep("S[^v]", names(res1))])]
        }                          
    }
    data.table::setattr(x, "sf", res1)
    return(attr(x, "sf"))
}
#' 
sf.numeric <- function(x, ..., n,
                       what=c("all", "S", "Sv", "H", "Hv")){
    if (missing(n)) sf.default(x)
    stopifnot(is.numeric(x) && is.numeric(n))
    stopifnot(length(x)==length(n))
    what <- match.arg(what)
    if(what == "all"){
        res1 <- data.table::data.table(
            "S"=km(x, n),
            "Sv"=kmv(x, n),
            "H"=na(x, n),
            "Hv"=nav(x, n))
        res1[, "sna" := exp(-H)]
        res1[, "hkm" := -log(S)]
        return(res1)
    }
    res1 <- switch(what,
                   "S"=km(x, n),
                   "Sv"=kmv(x, n),
                   "H"=na(x, n),
                   "Hv"=nav(x, n))
    return(res1)
}
###
### functions called above
###
km <- function(e, n) cumprod(1 - (e / n))
kmv <- function(e, n) cumprod(1 - (e / n))^2 * cumsum(e / (n * (n - e)))
na <- function(e, n) cumsum(e / n)
nav <- function(e, n) cumsum(e / (n^2))
