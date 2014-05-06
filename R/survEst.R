##' @name survEst
##' @aliases KM
##' @aliases Gw
##' @aliases NelA
##' @aliases NelAVar
##' @aliases survEst
##'
##' @title Estimates of survival functions based on n and e
##'
##' @param n Vector with corresponding no. at risk at time points
##' @param e Vector with no. events (e.g. deaths) at time points
##'
##' @note \code{survEst} is a generic name for the functions documented.
##'
##' @rdname survEst
##' @export
survEst <- function(n,e) {identity(c(n,e))}
##'
##' @rdname survEst
##' @include tne.R
##' @details \code{KM} gives the Kaplan-Meier (product-limit) estimator of survival function.
##' Value is \eqn{1} at \eqn{t=0} and thereafter given by:
##' \deqn{\prod_{t \leq t_i} (1-\frac{e_i}{n_i} )}{PROD 1-e/n}
##' @return \code{KM} returns \eqn{\hat{S}}{S}, the Kaplan-Meier estimator
##' \cr
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' t1 <- tne(s1)
##' KM(n=t1$n, e=t1$e)
##' ###
##' @export
KM <- function(n,e){
    if (!length(e)==length(n)) stop("Vectors must be of equal length")
    if (!is.numeric(e)|!is.numeric(n)) stop("Vectors must be numeric")
    cumprod(1-(e/n))
}
##' @rdname survEst
##' @include tne.R
##' @details \code{Gw} gives Greenwoods estimtor of the Kaplan-Meier (product-limit)
##' estimator of survival function is given by:
##'
##' \deqn{Var[\hat{S}(t)] = [\hat{S}(t)]^2 \sum_{t_i \leq t} \frac{e_i}{n_i (n_i-e_i)} }{
##' Var[S](t) = S(t)^2 SUM (e) / n(n-e) }
##' @return \code{Gw} returns \eqn{Var [\hat{S}(t)]},
##' the variance of the Kaplan-Meier estimator
##' \cr
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' t1 <- tne(s1)
##' Gw(n=t1$n, e=t1$e)
##' data(drug6mp, package="KMsurv")
##' s1 <- Surv(time=drug6mp$t2, event=drug6mp$relapse) # 6MP patients
##' t1 <- tne(s1)
##' Gw(n=t1$n, e=t1$e)
##' ###
##' @references Last example for \code{Gw} is from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Table 4.1A, pg 93.
##'
##' @export
Gw <- function(n,e){
        if (!length(e)==length(n)) stop("Vectors must be of equal length")
        if (!is.numeric(e)|!is.numeric(n)) stop("Vectors must be numeric")
        return(cumprod(1-(e/n))^2*  cumsum( e / (n*(n-e)) ) )
    }
##' @rdname survEst
##' @include tne.R
##' @details \code{NelA} gives the Nelson Aalen estimator of hazard function.
##' Value is \eqn{\hat{H}=0}{H=0} at \eqn{t=0} and thereafter given by:
##' \deqn{\hat{H}(t) = \sum_{t \leq t_i} \frac{e_i}{n_i}  }{ H(t) = SUM e/n}
##' @return \code{NelA} returns \eqn{\hat{H}(t)}, the Nelson-Aalen estimator
##' \cr
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' t1 <- tne(s1)
##' NelA(e=t1$e, n=t1$n)
##' ###
##' @export
NelA <- function(n,e){
    if (!length(e)==length(n)) stop("Vectors must be of equal length")
    if (!is.numeric(e)|!is.numeric(n)) stop("Vectors must be numeric")
    return(cumsum(e/n))
}
##' @rdname survEst
##' @export
##' @include tne.R
##' @details
##' \code{NelAVar} gives the variance of the Nelson-Aalen estimator,
##' given by: \deqn{ Var[\hat{H}(t))] = \sum_{t_i \leq t} \frac{e_i}{n_i^2} }{
##' SUM e/(n^2)}
##' @return \code{NelAV} returns \eqn{Var \hat{H}(t)}, the variance of
##' the Nelson-Aalen estimator
##' @examples
##' data(drug6mp, package="KMsurv")
##' s1 <- Surv(time=drug6mp$t2, event=drug6mp$relapse) # 6MP patients
##' t1 <- tne(s1)
##' NelAVar(n=t1$n, e=t1$e)
##' @references Example for \code{NelAVar} is from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Table 4.2, pg 94.
##' @export
NelAVar <- function(n,e){
    if (!length(e)==length(n)) stop("Vectors must be of equal length")
    if (!is.numeric(e)|!is.numeric(n)) stop("Vectors must be numeric")
    return(cumsum(e/(n^2)))
}
