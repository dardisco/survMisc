##' @name quantile
##' @export quantile.Surv
##' @aliases quantile.Surv
##' @method quantile Surv
##'
##' @include calcSurv.R
##' @include tne.Surv.R
##' @include survEst.R
##' @title Quantiles for \code{Surv} object
##' @param x A \code{Surv} object
##' @param ... Additional arguments
##' @param q Vector of quantiles (expressed as percentage)
##' @param alpha Significance level \eqn{\alpha}{alpha}
##'
##' @return A matrix with quantile, and upper and lower confidence intervals.
##' \cr
##' Intervals are calculated from \eqn{\sigma}{sigma} which is:
##' \deqn{\sigma (t) = \sqrt{ \frac{Var[\hat{S}(t)]}{\hat{S}^2(t)}} }{ ( Var[S(t)] / (S(t)^2) )^0.5}
##' The intervals given are:
##' \item{linear}{
##' \deqn{ \hat{S}(t) \pm Z_{1- \alpha} \sigma (t) \hat{S}(t)}{
##' S(t)+- Z(1-alpha) sigma(t) S(t)}
##' Where \eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival estimate.
##' \cr
##' }
##' \item{log transform}{
##' \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
##' [S(t)^(1/theta), S(t)^theta]}
##' Where \eqn{\theta}{theta} is:
##' \deqn{ e^{ \frac{Z_{1- \alpha} \sigma (t)}{ \log{\hat{S}(t)}}} }{
##' exp ( Z(1-alpha)sigma(t) / log(S(t)) )}
##' }
##' \item{arcsin-sqrt}{
##' Arcsine-square root transform.
##' \cr
##' \cr
##' Upper:
##' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}} - \frac{Z_{1- \alpha}\sigma(t)}{2} \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##' sin^2( max[0, arcsin S(t)^0.5 - Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##' Lower:
##' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}} + \frac{Z_{1- \alpha}\sigma(t)}{2} \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
##'  sin^2( min[pi/2, arcsin S(t)^0.5  + Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
##' }
##' @examples
##' data(bmt, package="KMsurv")
##' b1 <- bmt[bmt$group==1, ] # ALL patients
##' s1 <- Surv(time=b1$t2, event=b1$d3)
##' quantile(s1)
##'
quantile.Surv <- function(x, ..., q=c(25,50,75), alpha=0.05){
    if(!class(x)=="Surv") stop("Only applies to class 'Surv'")
    if(!attr(x,which="type")=="right") warning("Applies to right censored data")
###
    s2 <- x[order(x[,"time"]),]
    c1 <- calcSurv(x)[ ,c("t","SKM","SKMV")]
### express as percentage
    p1 <- q/100
### find corresponding time
    findT <- function(p) {
        ifelse(any(c1[,"SKM"]<=p),
               min(c1[which(c1[,"SKM"]<=p),"t"]),
               NA)
        }
    t1 <- sapply(p1,findT)
### convert alpha to z value
    z1 <- stats::qnorm(1-alpha/2)
### linear
    findLin <- function(p) {
        (c1[,"SKM"]-p) / sqrt((c1[,"SKMV"]))
            }
### log transform (negative log log)
    findLog <- function(p){
        ( ( log(-log(c1[,"SKM"])) - (log(-log(p))) )*c1[,"SKM"]*log(c1[,"SKM"]))/ sqrt(c1[,"SKMV"])
        }
### arcsin-square root transform
    findArc <- function(p){
        ( 2*(asin(sqrt(c1[,"SKM"])) - asin(sqrt(p)))*sqrt( c1[,"SKM"]*(1-c1[,"SKM"])) )/sqrt(c1[,"SKMV"])
    }
###
    Lin1 <- sapply(p1,findLin)
    Log1 <- sapply(p1,findLog)
    Arc1 <- sapply(p1,findArc)
    findMin <- function(mat){
        sapply(1:ncol(mat),
               function (i){
### if any values < z
                   ifelse(any(mat[,i]<=z1),
### get longest (largest) corresponding time
                          c1[min(which(mat[,i]<=z1)),"t"],
                          NA)
               }
               )
    }
    findMax <- function(mat){
        sapply(1:ncol(mat),
               function (i){
                   ifelse(any(mat[,i]<=-z1),
                          c1[max(which(mat[,i]<=-z1)),"t"],
                          NA)
               }
               )
    }
    res1 <- rbind(t1,
                findMin(Lin1),
                findMax(Lin1),
                findMin(Log1),
                findMax(Log1),
                findMin(Arc1),
                findMax(Arc1)
                )
    colnames(res1) <- q
    rownames(res1) <- c("quantile",
                      "lin.lower","lin.upper",
                      "log.lower","log.upper",
                      "a.s.lower","a.s.upper")
    return(res1)
}
