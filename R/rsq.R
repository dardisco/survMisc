##' @name rsq
##' @title R^2 measures for a a \code{coxph} or \code{survfit} model
##' @rdname rsq
##' @export rsq
##'
rsq <- function(x, ...){
    UseMethod("rsq")
}
##' @rdname rsq
##' @aliases rsq.coxph
##' @method rsq coxph
##' @S3method rsq coxph
##'
##' @param x A \code{survfit} or \code{coxph} object
##' @param ... Additional arguments (not implemented)
##' @param sigd Significant digits (for display). If \code{sigd=NULL}
##' will return the original numbers
##' @return A \code{list} with the following elements:
##'  \item{cod}{The coefficient of determination, which is
##'   \deqn{R^2=1-\exp(\frac{2}{n}L_0-L_1)}{
##'   R^2 = 1-exp((2/n).(L[0]-L[1]))}
##' where \eqn{L_0}{L[0]} and \eqn{L_1}{L[1]} are the log partial
##' likelihoods for the \emph{null} and \emph{full} models respectively
##' and \eqn{n}
##' is the number of observations in the data set.}
##' \item{mer}{The measure of explained randomness which is
##'   \deqn{R^2_{mer}=1-\exp(\frac{2}{m}L_0-L_1)}{
##'   R^2[mer] = 1-exp((2/m).(L[0]-L[1]))}
##' where \eqn{m} is the number of observed \emph{events}.}
##' \item{mev}{The measure of explained variation (similar to
##' that for linear regression) which is
##'   \deqn{R^2=\frac{R^2_{mer}}{R^2_{mer} + \frac{\pi}{6}(1-R^2_{mer})}}{
##'   R^2 = R^2[mer] / ( R^2[mer] + pi/6(1-R^2[mer]) )}
##' }
##'
##' @references Nagelkerke NJD, 1991.
##' A Note on a General Definition of the Coefficient of Determination.
##' \emph{Biometrika} \bold{78}(3):691--92.
##' \href{http://www.jstor.org/stable/2337038}{JSTOR}
##' @references
##' O'Quigley J, Xu R, Stare J, 2005.
##' Explained randomness in proportional hazards models.
##' \emph{Stat Med} \bold{24}(3):479--89.
##' \href{http://www.math.ucsd.edu/~rxu/igain2.pdf}{UCSD}
##' @references
##' Royston P, 2006.
##' Explained variation for survival models.
##' \emph{The Stata Journal} \bold{6}(1):83--96.
##' \href{http://www.stata-journal.com/sjpdf.html?articlenum=st0098}{The Stata Journal}
##'
rsq.coxph <- function(x, ..., sigd=2){
    stopifnot(inherits(x, "coxph"))
    l0 <- x$loglik[1]
    l1 <- x$loglik[2]
    n1 <- x$n
    ne1 <- x$nevent
    res1 <-  vector("list", length=3L)
    names(res1) <- c("nag", "mer", "mev")
    res1$cod <- 1 - exp((2 / n1) * (l0 - l1))
    res1$mer <- 1 - exp((2 / ne1) * (l0-l1))
    res1$mev <- res1$mer / (res1$mer + pi^2 / 6 * (1 - res1$mer))
    if (is.null(sigd)) return(res1)
    res1 <- lapply(res1, function(N)
                   as.numeric(formatC(signif(N, digits=sigd),
                                      digits=sigd, format="fg",
                                      flag="#")))
}
##' @rdname rsq
##' @aliases rsq.survfit
##' @method rsq survfit
##' @S3method rsq survfit
rsq.survfit <- function(x, ..., sigd=2){
    c1 <- deparse(x$call)
    c1 <- sub("survfit", "coxph", c1)
    c1 <- eval(parse(text=c1))
    rsq(c1, sigd=sigd)
}

