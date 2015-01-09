##' @name as.Surv
##' @rdname asSurv
##' @export as.Surv
##' @aliases as.Surv
##'
##' @title Convert \code{time} and \code{status} to a
##' right-censored survival object, ordered by time
##'
##' @param ti Vector of time points
##' @param st Vector of status (e.g. death) at time points
##' 
##' @return
##' An object of class \code{Surv}.
##'
##' @details
##' A traditional \code{Surv} object will only allow discrete events,
##' i.e. status \eqn{s \in N}{s = Natural number}.
##' Typically, \eqn{s \in \{0,1\}}{s = 0|1}.
##' \cr
##' There may be times when allowing non-binary values is of interest,
##' e.g. in constructing \emph{expected} survival. 
##' \cr
##' Caution is required when using this function with functions that
##' assume binary values for \code{status}.
##'
##' @seealso
##' \code{\link{dx}}
##' \cr
##' \code{?survival::Surv}
##'
##' @examples
##' c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
##' E <- predict(c1, type="expected")
##' as.Surv(pbc$time, E)
##' \dontrun{
##' summary(coxph(as.Surv(pbc$time, E) ~ log(pbc$bili)))
##' ### Warning:
##' ### In Surv(stime, sstat) : Invalid status value, converted to NA}
as.Surv <- function(ti, st){
    res1 <- cbind(ti, st)[order(ti), ]
    res1 <- res1
    colnames(res1) <- c("time", "status")
    class(res1) <- "Surv"
    attr(res1, "type") <- "right"
    return(res1)
}
