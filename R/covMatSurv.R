##' @name covMatSurv
##' @export
##' @title Covariance matrix for survival data
##' @description Gives variance-covariance matrix for comparing survival
##' data for two or more groups.
##' Inputs are vectors corresponding to observations at a set of discrete
##' time points for right censored data, except for \eqn{n1},
##' the no. at risk by predictor.
##' This should be specified as a vector for one group,
##' otherwise as a matrix with each column corresponding to a group.
##' @param t time
##' @param n no. at risk
##' @param e no. events
##' @param n1 no. at risk (by predictor).
##' \cr
##' If 2 groups, should be given as a vector n1 with no. at risk for n1.
##' \cr
##' If more than 2 groups, a matrix with a column for each group.
##' @return An array. First two dimensions = no. groups.
##' \cr
##' Third dimension is no. observations (time points).
##' \cr \cr
##' Where there are two groups, the resulting sparse square matrix
##' at time \eqn{i} has diagonal elements:
##'  \deqn{v_i = - \frac{n0_i n1_i e_i (n_i-e_i)}{n_i^2(n_i-1)}}{
##' n0(i).n1(i).e(i).(n(i)-e(i)) / n(i)^2.(n(i) -1)}
##' where \eqn{n1} is the no. at risk in group 1.
##' \cr \cr
##' For more than two groups, the resulting square matrix has diagonal elements:
##' \deqn{ v_{kki} = \frac{n_{ki}(n_i-n_{ki})e_i(n_i-e_i)}{n_i^2(n_i-1)}}{
##'  v[k,k,i] = n[k](i).[n(i) - n[k](i)].e(i).[n(i) -e(i)] / n(i)^2.(n(i)-1)}
##' and off diagonal elements:
##' \deqn{ v_{kli} = \frac{ -n_{ki}n_{li} e_i(n_i-e_i)}{n_i^2(n_i-1)}}{
##'  v[k,l,i] = - n[k](i).n[l](i).e(i).[n(i)-e(i)] / n(i)^2.[n(i)-1]}
##' @seealso Called by \code{\link{comp}}
##' @examples
##' \dontrun{data(tneKidney)
##' covMatSurv(t=tneKidney$t,n=tneKidney$n,e=tneKidney$e,n1=tneKidney$n_1)
##' }
covMatSurv <- function(t,n,e,n1){
    if(!isTRUE( all(length(t)==length(n),
                    length(t)==length(e),
                    length(t==nrow(as.matrix(n1)))))) stop ("All vectors must be of equal length")
    if(!isTRUE(all(vapply(c(t,n,e,n1),FUN=is.numeric,FUN.VALUE=TRUE)==TRUE))) stop("All vectors must be numeric")
### no. of groups
    g1 <- ifelse(is.null(dim(n1)),1,ncol(n1))
### if 2 groups only
    if (g1==1){
        var1 <- (n1/n)*(1-(n1/n))*((n-e)/(n-1))*e
        return(var1)
        }
### hold results
    a1 <- array(data=0, dim=c(g1,g1,length(t)))
### diagonal elements
    for (i in 1:g1){
        a1[i,i, ] <- ( n1[,i]*(n-n1[,i])*e*(n-e) ) / ( (n^2)*(n-1) )
    }
### off-diagonal elements
    for (j in 1:g1){
        for (k in 1:g1){
           if(j==k) next
           a1[j,k, ] <- -( n1[,j]*n1[,k]*e*(n-e) ) / ( (n^2)*(n-1) )
           }
        }
    dimnames(a1) <- list(1:g1,1:g1,t)
    return(a1)
}
