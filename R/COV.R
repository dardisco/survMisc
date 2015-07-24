#' @name COV
#' @title \bold{Cov}ariance matrix for survival data
#'
#' @include tn.R
#' 
#' @rdname COV
#' @export
#'
COV <- function(x, ...) UseMethod("COV")
#'
#'
#' @param x
#'  For the default method, a \code{numeric} vector of
#'  \emph{number of events}, \eqn{e_t}{e[t]}.
#'  These are assumed to be ordered by discrete times.
#'  \cr
#'  A method is available for objects of \code{class} \code{tne}.
#' @param ... Additional arguments (not implemented).
#'  \cr
#'  The following arguments apply only to the default method.
#' @param n \bold{N}umber at risk (total).
#' @param ncg \bold{N}umber at risk, per \bold{c}ovariate \bold{g}roup.
#'  \cr
#' If there are \eqn{2} groups, this can be given as a \code{vector} with
#' the number at risk for group \eqn{1}.
#' \cr
#' If there are \eqn{\geq 2}{>=2} groups, it is
#' a \code{matrix} with one column for each group.
#'
#' @details Gives variance-covariance matrix for comparing survival
#' data for two or more groups.
#' \cr
#' Inputs are vectors corresponding to observations at a set of discrete
#' time points for right censored data, except for \eqn{n1},
#' the no. at risk by predictor.
#' \cr
#' This should be specified as a vector for one group,
#' otherwise as a matrix with each column corresponding to a group.
#'
#' @return An \code{array}.
#' \cr
#' The first two dimensions = the number of covariate groups \eqn{K},
#' \eqn{k = 1, 2, \ldots K}.
#' This is the square matrix below.
#' \cr
#' The third dimension is the number of observations
#' (discrete time points).
#' \cr \cr
#' To calculate this, we use \code{x} (= \eqn{e_t}{e[t]} below) and
#' \eqn{n_1}{n1}, the number at risk in covariate group \eqn{1}.
#' \cr
#' Where there are \eqn{2} groups, the resulting sparse square matrix
#' (i.e. the non-diagonal elements are \eqn{0})
#' at time \eqn{t} has diagonal elements:
#'  \deqn{cov_t = - \frac{n_{0t} n_{1t} e_t (n_t - e_t)}{n_t^2(n_t-1)}}{
#'        cov[t] = - n0[t] * n1[t] * e[t] * (n[t] - e[t]) /
#'                  (n[t]^2 * (n[t] - 1))}
#' For \eqn{\geq 2}{>=2} groups, the resulting square matrix
#' has diagonal elements given by:
#'  \deqn{cov_{kkt} = \frac{n_{kt}(n_t - n_{kt}) e_t(n_t - e_t)}{
#'                          n_t^2(n_t - 1)}}{
#'    cov[k, k, t] = n[k, t] * (n[t] - n[k, t]) * e[t] * (n[t] - e[t]) /
#'                   (n[t]^2 * (n[t] - 1))}
#' The off diagonal elements are:
#' \deqn{cov_{klt} = \frac{-n_{kt} n_{lt} e_t (n_t-e_t) }{
#'                         n_t^2(n_t-1)}}{
#'       cov[k, l, t] = - n[k, t] * n[l, t] * e[t] * (n[t] - e[t]) /
#'                      n[t]^2 * (n[t] - 1)}
#'
#' @seealso Called by \code{\link{comp}}
#'
#' @keywords survival
#'
#' @rdname COV
#' @aliases COV.tn
#' @method COV tn
#' @export
#'
#' @examples
#' data(kidney, package="KMsurv")
#' k1 <- with(kidney,
#'             tne(Surv(time=time, event=delta) ~ type)
#' with(k1, COV(x=e, n=n, n1=n_1))
#' COV(k1)
#'
COV.tn <- function(x, ...){
    stopifnot(attr(x, "byWhat")=="time")
    ## no. of groups
    g1 <- attr(x, "ncg")
    if (g1 <= 1) error("Only valid if more than one covariate group")
    ## if 2 groups only
    if (g1==2){
      n1 <- ifelse(attr(x, "short.Names"),
                   "n_1",
                   paste0)
      n1 <- names(x)
      data.table::setnames(x, c(names(x)[1:3],
                                "n_1",
                                names(x)[5:length(names(x))]))
      data.table::setattr(x, "cov",
                          x[,
                            (n_1 / n) * (1 - (n_1 / n)) *
                            ((n - e) / (n - 1)) * e]
    }
    ## more than 2 groups?
    if(g1 > 2){
    ### hold results
        res1 <- array(data=0,  dim=c(g1, g1, nrow(x)))
### diagonal elements
        for (i in seq.int(g1)){
            res1[i, i, ] <- (n1[, i] * (n - n1[, i]) * e * (n - e)) / ((n^2) * (n - 1))
        }
### off-diagonal elements
        for (j in seq.int(g1)){
            for (k in seq.int(g1)){
                if (j==k) next
                res1[j, k, ] <- -( n1[, j]*n1[, k]*e*(n-e) ) / ( (n^2)*(n-1) )
            }
        }
        dimnames(res1) <- list(1:g1, 1:g1, t)
    }
    return(res1)
      }
    
COV.stratTn <- function(x, ...){
    return(lapply(x, FUN=coV))
}
###----------------------------------------
## covMatSurv <- function(t, n, e, n1){
##     stopifnot(all(sapply(list(t, n, e, n1), is.numeric)))
##     ## ensure all same length
##     stopifnot(diff(range(sapply(list(t, n, e), length))) < .Machine$double.eps)
## ### no. of groups
##     g1 <- ifelse(
##         is.null(dim(n1)) | ncol(n1)==1,
##         1,
##         ncol(n1))
## ### if 2 groups only
##     if (g1==1){
##         var1 <- (n1 / n) * (1 - (n1 / n)) * ((n - e) / (n - 1)) * e
##         return(var1)
##     }
## ### hold results
##     a1 <- array(data=0,  dim=c(g1, g1, length(t)))
## ### diagonal elements
##     for (i in seq_len(g1)){
##         a1[i, i, ] <- (n1[, i] * (n - n1[, i]) * e * (n - e)) / ((n^2) * (n - 1))
##     }
## ### off-diagonal elements
##     for (j in seq_len(g1)){
##         for (k in 1:g1){
##             if (j==k) next
##             a1[j, k, ] <- -( n1[, j]*n1[, k]*e*(n-e) ) / ( (n^2)*(n-1) )
##         }
##     }
##     dimnames(a1) <- list(1:g1, 1:g1, t)
##     return(a1)
## }
