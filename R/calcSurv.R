##' @name calcSurv
##' @export
##' @include tne.Surv.R
##' @include survEst.R
##' @title Estimators of survival & hazard functions with variance for a \code{Surv} object describing right censored data
##' @param s A \code{Surv} object
##' @return A data frame. \cr \cr
##' Rows - one for each time where an event occurs.
##' \cr \cr
##' Columns are:
##'  \item{t}{time}
##'  \item{n}{no. at risk}
##'  \item{e}{no. events}
##'  \item{KM}{Survival estimate by Kaplan-Meier (Product-Limit) estimator}
##'  \item{KMV}{Variance of Kaplan-Meier estimate (Greenwoods formula)}
##'  \item{SNel}{
##'   Survival estimate from Nelson-Aalen estimator:
##'    \eqn{\hat{S}=e^{\hat{H}} }{S=e^H}
##'  }
##'  \item{HNel}{Nelson-Aalen estimate of hazard function}
##'  \item{HNelV}{Variance of Nelson-Aalen estimate}
##'  \item{HKM}{
##'   Hazard estimate from Kaplan-Meier estimator:
##'    \eqn{\hat{H}=-log{\hat{S} }}{H = -log(S)}
##'   }
##' @seealso \code{\link{tne.Surv}}
##' @seealso \code{\link{KM}}
##' @seealso \code{\link{Gw}}
##' @seealso \code{\link{NelA}}
##' @seealso \code{\link{NelAVar}}
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' calcSurv(s1)
calcSurv <- function(s){
    if(!class(s)=="Surv") stop("Only applies to class 'Surv'")
    if(!attr(s,which="type")=="right") warning("Applies to right censored data")
### sort by time (makes debugging easier)
    s2 <- s[order(s[,"time"]), ]
### hold results
    m1 <- data.frame(matrix(nrow=sum(s[,"status"]),ncol=9))
### if last observation censored then generate final row
    if(!s2[,"status"][s2[,"time"]==max(s2[,"time"])]){
        s3 <- s2[which(s2[,"time"]==max(s2[,"time"]))]
        tne1 <- tne.Surv(s2,onlyEvents=FALSE)
        m1[nrow(m1),1] <- tne1[nrow(tne1),1]
        m1[nrow(m1),2] <- tne1[nrow(tne1),2]
        m1[nrow(m1),3] <- tne1[nrow(tne1),3]
### time, no, event
        m1[1:(nrow(m1)-1),1:3] <- tne.Surv(s)
    } else {
        m1[,1:3] <- tne.Surv(s2)
        }
### names to fill in
    colnames(m1) <- c("t","n","e",
                      "SKM","SKMV","SNelA",
                      "HNel","HNelV","HNelV")
### Kaplan-Meier
    m1$SKM <- KM(e=m1$e,n=m1$n)
### Greenwoods formula for variance
    m1$SKMV <- Gw(e=m1$e,n=m1$n)
    m1$SNelA <- exp(-(NelA(e=m1$e,n=m1$n)))
### Nelson-Aalen
    m1$HNel <- NelA(e=m1$e,n=m1$n)
### variance
    m1$HNelV <- NelAVar(e=m1$e,n=m1$n)
    return(m1)
}
