#' @name covMatSurv
#' @export
#' @title Covariance matrix for survival data
#' 
#' @param t time.
#' @param n number at risk (total).
#' @param e number of events (total).
#' @param ncc number at risk (by predictor).
#' \cr
#' If \eqn{2} groups, should be given as a \code{vector} with
#' the number at risk for group \eqn{1}.
#' \cr
#' If \eqn{\geq 2}{>=2} groups, a \code{matrix} with one column for each group.
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
#' @return An \code{array}. The first two dimensions = number of groups.
#' This is the square matrix below.
#' \cr
#' The third dimension is the number of observations (time points).
#' \cr \cr
#' Where there are \eqn{2} groups, the resulting sparse square matrix
#' (i.e. the non-diagonal elements are \eqn{0})
#' at time \eqn{i} has diagonal elements:
#'  \deqn{v_i = - \frac{n_{0i} n_{1i} e_i (n_i-e_i)}{n_i^2(n_i-1)}}{
#'        v(i) = - n0(i).n1(i).e(i).(n(i)-e(i)) / n(i)^2.(n(i) -1)}
#' where \eqn{n_1}{n1} is the number at risk in group \eqn{1}.
#' \cr \cr
#' For \eqn{\geq 2}{>=2} groups, the resulting square matrix has diagonal elements:
#' \deqn{ v_{kki} = \frac{n_{ki}(n_i - n_{ki})e_i(n_i - e_i)}{n_i^2(n_i - 1)}}{
#'        v[k, k, i] = n[k, i] * n[i] - n[k, i] * e[i] * (n[i] - e[i]) / n[i]^2 * (n[i] - 1)}
#' and off diagonal elements:
#' \deqn{ v_{kli} = \frac{-n_{ki}n_{li} e_i(n_i-e_i)}{n_i^2(n_i-1)}}{
#'        v[k, l, i] = - n[k, i]* n[l, i] * e[i] * (n[i] - e[i]) / n[i]^2 * (n[i] - 1)}
#' 
#' @seealso Called by \code{\link{comp}}
#'
#' @keywords survival
#'
#' @examples
#' data(package="KMsurv")
#' data(kidney, package="KMsurv")
#' k1 <- with(kidney,
#'            tne(Surv(time=time, event=delta) ~ type, shortNames=TRUE))
#' with(k1, covMatSurv(t, n=n, e, n1=n_1))
covMatSurv <- function(t, n, e, n1){
    stopifnot(all(sapply(list(t, n, e, n1), is.numeric)))
    ## ensure all same length
    stopifnot(diff(range(sapply(list(t, n, e), length))) < .Machine$double.eps)
### no. of groups
    g1 <- ifelse(
        is.null(dim(n1)) | ncol(n1)==1,
        1,
        ncol(n1))
### if 2 groups only
    if (g1==1){
        var1 <- (n1 / n) * (1 - (n1 / n)) * ((n - e) / (n - 1)) * e
        return(var1)
    }
### hold results
    a1 <- array(data=0,  dim=c(g1, g1, length(t)))
### diagonal elements
    for (i in seq_len(g1)){
        a1[i, i, ] <- (n1[, i] * (n - n1[, i]) * e * (n - e)) / ((n^2) * (n - 1))
    }
### off-diagonal elements
    for (j in seq_len(g1)){
        for (k in 1:g1){
            if (j==k) next
            a1[j, k, ] <- -( n1[, j]*n1[, k]*e*(n-e) ) / ( (n^2)*(n-1) )
        }
    }
    dimnames(a1) <- list(1:g1, 1:g1, t)
    return(a1)
}
