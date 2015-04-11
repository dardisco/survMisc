##' @name ic
##' @rdname ic
##' @title Information criterion
##' @aliases AIC
##' @aliases Akaike information criterion
##' @aliases Bayesian information criterion
##' @aliases AICc
##' @aliases BIC
##' 
##' @param object An object of class \code{coxph}
##' @param ... Not implemented
##' @param k The weight of the equivalent degrees of
##' freedom (edf) of the AIC formula
##' @description \bold{I}nformation \bold{C}riterion for a fitted model.
##' 
##' @details Given a set of candidate models for the same data,
##' the preferred model is the one with the minimum IC value.
##' \cr
##' The Akaike information criterion, AIC, is given by
##' \deqn{AIC = k.edf -2 \ln L}{
##'  AIC = k*edf - 2*log L}
##' Where \eqn{edf} is the
##' equivalent degrees of freedom
##' (i.e., equivalent to the number of free parameters in the model) and
##' \eqn{L} is the model likelihood.
##' \cr
##' \eqn{k} is a constant, which is \eqn{=2} for the traditional AIC.
##' \cr \cr
##' AIC corrected for finite sample size \eqn{n}, AICc, is
##' \deqn{AICc = AIC + \frac{k.edf(edf+1)}{n-edf-1}}{
##'  AICc = AIC + k*edf*(edf+1) / (n-edf-1)}
##' where \eqn{n} is the sample size. Thus there is a greater penalty
##' for more parameters.
##' \cr \cr
##' The Bayesian information criterion is
##' \deqn{BIC = \ln n.edf -2 \ln L}{
##'  BIC = log(n)*edf - 2*log L}
##' This  penalises models with more parameters to a greater extent.
##' 
##' @note For survival models the \strong{effective \eqn{n}}
##' is the number of events rather
##' than the number of observations. This is used in computing the criteria
##' above.
##' 
##' @return A named vector with \describe{
##'   \item{edf}{the equivalent degrees of freedom for the fitted model fit}
##'   \item{ IC}{the information criterion, either AIC, AICc or BIC}
##' }
NULL
##' @rdname ic
##' @export
BIC <- function(object, ...){
    UseMethod("BIC")
}
##' @rdname ic
##' @method BIC coxph
##' @export
BIC.coxph <- function(object, ...){
    n1 <- object$nevent
    a1 <- stats::extractAIC(object, k=log(n1))
    names(a1) <- c("edf", "BIC")
    return(a1)
}
##' @rdname ic
##' @export
AIC <- function(object, ..., k=2){
    UseMethod("AIC")
}
##' @rdname ic
##' @method AIC coxph
##' @export
AIC.coxph <- function(object, ..., k=2){
    a1 <- stats::extractAIC(object, k)
    names(a1) <- c("edf", "AIC")
    return(a1)
}
##' @rdname ic
##' @export
AICc <- function(object, ...){
    UseMethod("AICc")
}
##' @rdname ic
##' @method AICc coxph
##' @export
AICc.coxph <- function(object, ..., k=2){
    n1 <- object$nevent
### correction
    corr1 <- 2*k*(k+1) / (n1-k-1)
    a1 <- stats::extractAIC(object, k)
    a1[2] <- a1[2] + corr1
    names(a1) <- c("edf", "AICc")
    return(a1)
}
