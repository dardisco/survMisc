##' @name calcSurv
##' @export
##' @include tne.R
##' @include survEst.R
##' @title Calculate estimators of survival & hazard functions
##' with variance for a \code{Surv} object describing right censored data
##' @param s A \code{Surv} object
##' @param round No. digits to which to round (for display)
##' @return A data frame. \cr \cr
##' Rows - one for each time where an event occurs.
##' \cr \cr
##' Columns are:
##'  \item{t}{time}
##'  \item{n}{no. at risk}
##'  \item{e}{no. events}
##'  \item{KM}{Survival estimate by Kaplan-Meier (Product-Limit) estimator}
##'  \item{KMV}{Variance of Kaplan-Meier estimate (Greenwoods formula)}
##'  \item{SNelA}{
##'   Survival estimate from Nelson-Aalen estimator:
##'    \eqn{\hat{S}=e^{\hat{H}} }{S=e^H}
##'  }
##'  \item{HNel}{Nelson-Aalen estimate of hazard function}
##'  \item{HNelV}{Variance of Nelson-Aalen estimate}
##'  \item{HKM}{
##'   Hazard estimate from Kaplan-Meier estimator:
##'    \eqn{\hat{H}=-log{\hat{S} }}{H = -log(S)}
##'   }
##' @seealso \code{\link{tne}}
##' @seealso \code{\link{KM}}
##' @seealso \code{\link{Gw}}
##' @seealso \code{\link{NelA}}
##' @seealso \code{\link{NelAVar}}
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' calcSurv(s1)
##' data.frame(calcSurv(s1, round=4)[1:4],
##'           round(sqrt(calcSurv(s1)[5]),4),
##'           calcSurv(s1, round=4)[7],
##'           round(sqrt(calcSurv(s1)[8]),4)
##'           )
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==2, ] # AML patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' calcSurv(s1)
##' @references First example is from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Table 4.3, pg 97.
calcSurv <- function(s, round=7){
    if(!class(s)=="Surv") stop("Only applies to class 'Surv'")
    if(!attr(s,which="type")=="right") warning("Applies to right censored data")
### sort by time (makes debugging easier)
    s <- s[order(s[, "time"]), ]
    m1 <- tne(s)
### Kaplan-Meier
    m1$SKM <- KM(e=m1$e, n=m1$n)
### Greenwoods formula for variance
    m1$SKMV <- Gw(e=m1$e, n=m1$n)
### Nelson-Aalen (hazard)
    m1$HNel <- NelA(e=m1$e, n=m1$n)
### get survival survival from NelA estimate of hazard
    m1$SNelA <- exp(-m1$HNel)
### variance
    m1$HNelV <- NelAVar(e=m1$e, n=m1$n)
### get hazard from KM estimate of survival
    m1$HKM <- -log(m1$SKM)
    m1 <- round(m1, round)
    return(m1)
}
