#' @name sf
#' @export sf
#'
#' @title Estimates of \bold{s}urvival (or hazard) \bold{f}unction
#' based on \eqn{n} and \eqn{e}
#'
#' @include tne.R
#'
#' @param x
#' @param n \bold{N}umber at risk per time point (a \code{vector})
#' @param e Number of \bold{e}vents per time point (a \code{vector})
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
#' data(drug6mp, package="KMsurv")
#' s1 <- with(drug6mp, Surv(time=t2, event=relapse))
#' t1 <- tne(s1)
#' sf(n=t1$n, e=t1$e, what="sv")
#' ## K&M 2nd ed. Table 4.2, pg 94.
#' data(bmt, package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' t2 <- tne(Surv(time=b1$t2, event=b1$d3))
#' with(t2, sf(n=n, e=e, what="hv"))
#' ## K&M 2nd ed. Table 4.3, pg 97.
#' sf(n=t2$n, e=t2$e, what="all")
#'
sf <- function(x, ...) UseMethod("sf")

t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, events.Only=FALSE, by.What="time")

sf.default <- function(x, ...){
    stopifnot(all(x %in% c(0,1)))
    t1 <- tne(x)
    return(sf.tne(t1))
}
#'
#'
t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, events.Only=FALSE, by.What="time")
debugonce("sf.tne")
sf(t1)
sf.tne <- function(x, ...,
                   what=c("S", "H")){
    what <- match.arg(what)
    ## functions to use
    if(what=="S"){
        fun1 <- km
        fun1v <- kmv
    } else {
        fun1 <- na
        fun1v <- nav
    }
    if(attr(x, "by.What")=="status"){
        ## name of variance
        nv1 <- paste0(what, "v")
        x[, (what) := fun1(ncg, status), by=cg]
        x[, (nv1) := fun1v(ncg, status), by=cg]
    }
    if(attr(x, "by.What")=="time"){
        n1 <- r1 <- e1 <- grep("n_", names(x), value=TRUE)
        substr(r1, 1L, 2L) <- paste0(what, "_")
        substr(e1, 1L, 2L) <- "e_"
        r1v <- gsub("n_", paste0(what, "v_"), n1)
        for(i in seq.int(attr(x, "ncg"))){
            x[, (r1)[i] := fun1(
                        eval(as.name(n1[i])),
                        eval(as.name(e1[i])))]
            x[, (r1v)[i] := fun1v(
                         eval(as.name(n1[i])),
                         eval(as.name(e1[i])))]
        }
        ## get names for new column order
        tne1 <- c("t", "n", "e")
        ne1 <- c(n1, e1)
        new1 <- c(r1, r1v)
        old1 <- names(object)[!names(object) %in%
                              c(tne1, ne1, new1)]
        ## co1 = column order
        co1 <- c(ne1, new1, old1)
        ## new column order
        co1 <- c(tne1, as.vector(t(matrix(co1, nrow=length(n1)))))
        data.table::setcolorder(object, co1)
    }
    data.table::setattr(x, "sf", what)
    return(x)
}
#'
sf.numeric <- function(x=NULL, ...,
                       n, e,
                       what=c("all", "s", "sv", "h", "hv")){
    res1 <- s <- sv <- h <- hv <- sna <- hkm <- NULL
    stopifnot(is.numeric(e) && is.numeric(n))
    stopifnot(length(e)==length(n))
    what <- match.arg(what)
    if(what == "all"){
        res1 <- data.table::data.table(s = km(n, e))
        res1[, "sv" := kmv(n, e)]
        res1[, "h" := na(n, e)]
        res1[, "hv" := nav(n, e)]
        res1[, "sna" := exp(-h)]
        res1[, "hkm" := -log(s)]
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
kmv <- function(n, e) cumprod(1 - (e / n))^2 * cumsum(e / (n * (n - e)))
na <- function(n, e) cumsum(e / n)
nav <- function(n, e) cumsum(e / (n^2))
