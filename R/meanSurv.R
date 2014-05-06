#' @name mean.Surv
#' @method mean Surv
#' @S3method mean Surv
#' @export
#' @include calcSurv.R
#' @title Mean for \code{Surv} object
#' @param x A \code{Surv} object
#' @param alpha Significance level \eqn{\alpha}{alpha}
#' @param method If the last observation is censored at time \eqn{t_k}{t_k}, one of the following values for \eqn{\hat{S}}{S}, the Kaplan-Meier estimate of survival time from then until \code{tMax} is used:
#' \describe{
#' \item{Effron}{ \eqn{\hat{S}=0}{S=0}}
#' \item{Gill}{ \eqn{\hat{S}=\hat{S}(t_k)}{S = S(t_k)} i.e. \eqn{\hat{S}}{S} is equal to the
#' last recorded value of \eqn{\hat{S}}{S}}.
#' \item{Brown}{ \eqn{\hat{S}= e^{\frac{t_i}{t_k} \log{\hat{S}(t_k)}} }{
#'  e^[ (t_i/t_k) log S(t_k) ]}
#' for \eqn{ t_k \le t_i \le \code{tMax}}{
#'  t(k) <= t(i) <= tMax}
#' }}
#'
#' @param tMax If the last observation is censored at time \eqn{t_k}{t_k}, an estimate of
#' \eqn{\hat{S}}{S} will be generated from \eqn{t_k}{t_k} to \code{tMax}.
#' \cr
#' If \code{tMax=NULL} a value of \eqn{2 \times t_{max}}{2*tMax}, twice the longest time recorded, is used.
#' @param by Increments (units of time) between \eqn{t_k}{t_k} and \code{tMax}
#' @param dfm If \code{TRUE}, will return the dataframe used to calculate the statistics
#' @param ... Additional arguments
#' @return A list with the following elements:
#' \item{mean}{Mean of the \code{Surv} object}
#' \item{variance}{The variance}
#' \item{CI}{The confidence level (from \eqn{alpha}{alpha} above)}
#' \item{upper}{Upper value for confidence interval}
#' \item{lower}{Lower value for the confidence interval}
#' If the last observation is censored at time \eqn{t_k}{t_k},
#' two values are returned, one calculated up to \eqn{t_k}{t_k}, the other to \code{tMax}.
#' @examples
#' data(bmt, package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' s1 <- Surv(time=b1$t2, event=b1$d3)
#' mean(s1)
#' mean(Surv(time=c(6, 14, 21, 44, 62), event=c(1, 1, 0, 1, 1)))
#' mean(Surv(time=c(6, 14, 21, 44, 62), event=c(1, 1, 0, 1, 0)))
#'
mean.Surv <- function(x, alpha=0.05,
                      method=c("Efron", "Gill", "Brown"),
                      tMax=NULL,
                      by=1,
                      dfm=FALSE, ...){
    if(!class(x)=="Surv") stop("Only applies to class 'Surv'")
    if(!attr(x,which="type")=="right") warning("Applies to right censored data")
    method <- match.arg(method)
### makes debugging easier
    s2 <- x[order(x[, "time"]), ]
### longest time
    t1 <- max(s2[, 1])
### longest observed time with an event
    t2 <- max(s2[, 1][s2[, 2]==1])
### calcSurv returns data.table by default
    c1 <- calcSurv(x)
    set(c1, j=5:9, value=NULL)
    class(c1) <- "data.frame"
### first row is t=0
    c1 <- rbind(c(0, NA, 0, 1), c1)
### no. at risk at t=0
    c1[1,"n"] <- c1[2,"n"]
###
### if last obseration censored add tail
    if (!t1==t2) {
### survival estimate for last observation
        SKMmin <- min(c1[,"SKM"])
### function for tail by method 'Brown' below
        BrownTail <- function(t) exp( (t/t1)*log(SKMmin) )
### extend to tMax
        if (is.null(tMax)) tMax <- t1*2
        t3 <- seq(from=t1,to=tMax,by=by)
### no. at risk at last observation
        nR <- sum(x[,"time"]==t1)
        tail1 <- switch(method,
                        Efron = matrix(c(t1, nR, 0, 0), nrow=1),
                        Gill = rbind(c(t1, nR, 0, SKMmin), c(tMax, nR, 0, SKMmin)),
                        Brown = cbind(t3,
                        rep(nR,length(t3)),
                        rep(0,length(t3)),
                        sapply(t3,BrownTail)),
                        )
        colnames(tail1) <- colnames(c1)
        c1 <- rbind(c1,tail1)
    }
    tMax <- max(c1[,"t"])
    rownames(c1) <- NULL
### area under curve = mean
    AUC <- function(i){
        (c1[i,"t"]-c1[(i-1),"t"])*c1[(i-1),"SKM"]
        }
### use tail (if censored)
    AUC1 <- sapply((nrow(c1)):2, AUC)
    auc1 <- cumsum(AUC1)
### exclude tail (if censored)
    n1 <- sum(c1[,"t"]<=t2)
    AUC1 <- sapply(n1:2, AUC)
    auc2 <- cumsum(AUC1)
### if largest observation censored,
### return 2x possible values of mean
   if (t1==t2){
       mean1 <- utils::tail(auc1,1)
       names(mean1) <- paste("t",t2,sep="")
       } else {
           mean1 <- c(utils::tail(auc1,1),
                      utils::tail(auc2,1))
           names(mean1) <- c(paste("t",tMax,sep=""),(paste("t",t2,sep="")))
           }
### variance using tail
    varf1 <- function(i){
      rev(auc1)[i]^2 * c1[i,"e"]/(c1[i,"n"]*(c1[i,"n"]-c1[i,"e"]))
      }
### variance excluding tail
    varf2 <- function(i){
      rev(auc2)[i]^2 * c1[i,"e"]/(c1[i,"n"]*(c1[i,"n"]-c1[i,"e"]))
      }
### if largest observation censored,
### return 2x possible values of variance
 if (t1==t2){
       var1 <- sum(sapply(1:(nrow(c1)-1),varf1))
       names(var1) <- paste("t",t2,sep="")
       } else {
           var1 <- c (
               sum(sapply(1:(nrow(c1)-1),varf1)),
               sum(sapply(1:(n1-1),varf2))
               )
           names(var1) <- c(paste("t",tMax,sep=""),
                            paste("t",t2,sep=""))
       }
### no. times with non-censored observations
    obs1 <- sum(c1$e>=1)
### Anderson's correction
    varC <-(var1*obs1) / (obs1 -1)
### change alpha to Z value
    z1 <- stats::qnorm(1-alpha/2)
### upper + lower
    u1 <- rbind(mean1 + z1*sqrt(var1),
                mean1 + z1*sqrt(varC))
    l1 <- rbind(mean1 - z1*sqrt(var1),
                mean1 - z1*sqrt(varC))
    rownames(u1) <- rownames(l1) <- c("unCorr", "Corr")
###
    CI <- paste(100*(1-alpha),"% CI",sep="")
    res <- list(
        "mean" = mean1,
        "variance" = var1,
        "CI" = CI,
        "upper" = u1,
        "lower" = l1)
    if (dfm) res$dfm <- c1
    return(res)
}
