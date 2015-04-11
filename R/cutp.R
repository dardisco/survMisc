##' @name cutp
##' @title Cutpoint for a continuous variable in a \code{coxph} or \code{survfit} model
##' @description Determine the optimal cutpoint for a continuous variable
##' in a \code{coxph} or \code{survfit} model
##' @rdname cutp
##' @export
##'
cutp <- function(x, ...){
    UseMethod("cutp")
    }
##' @rdname cutp
##' @aliases cutp.coxph
##' @method cutp coxph
##' @export
##'
##' @include tne.R
##' 
##' @param x A \code{survfit} or \code{coxph} object
##' @param ... Additional arguments. Passed to \code{graphics::plot}.
##' @param var Variable to test. Must be continuous (i.e. \eqn{>2} unique values)
##' @param plot If \code{plot=TRUE} will plot cut points against the test statistic \eqn{Q}.
##' @return A \code{data.frame} with columns:
##'  \item{cp}{The \bold{c}ut \bold{p}oint. The optimum value at which to divide
##' the groups into those \eqn{\geq}{>=} the cutpoint and those below.}
##'  \item{Q}{The test statistic}
##'  \item{p}{p-value}
##' If \code{plot=TRUE} a plot of cut points against values of the
##' log-rank test statistic \eqn{LR}.
##' 
##' @details
##' The statistic is based on the score test from the Cox model.
##' For the cut point \eqn{\mu}{mu}, of a predictor \eqn{K}, the data is split
##' into two groups, those \eqn{\geq \mu}{>= mu} and
##' those \eqn{< \mu}{< mu}.
##' \cr \cr
##' The log-rank statistic \eqn{LR} is calculated for each unique element
##' \eqn{k} in \eqn{K}:
##' \deqn{LR_k = \sum_{i=1}^D ( e_i^+ - n_i^+ \frac{e_i}{n_i} )}{
##'       LR[k] = sum ( e1[i] - n1[i].e[i]/n[i] ) }
##' Where \eqn{e_i^+}{e1[i]} and \eqn{n_i^+}{n1[i]} refer to the number of events
##' and number at risk in those above the cutpoint, respectively.
##' \cr
##' The sum is taken to across distinct times with observed events, to \eqn{D},
##' the largest of these.
##' \cr
##' It is normalized (standardized), in the case of censoring,
##' by finding \eqn{\sigma^2}{s^2} which is:
##' \deqn{ \sigma^2 = \frac{1}{D-1} \sum_i^D ( 1 - \sum_{j=1}^i \frac{1}{D+1-j} )^2 }{
##'        s^2 = 1/(D-1) SUM[i to D] { 1 - SUM[j to i] (1/(D-j+1))}^2 }
##' The test statistic is then
##' \deqn{Q = \frac{\max |LR_k|}{\sigma \sqrt{D-1}} }{
##'       Q = [ max |LR[k]| ] / [ s.(D-1)^0.5 ] }
##' Under the null hypothesis that the chosen cut-point does \emph{not} predict survival,
##' the distribution of \eqn{Q} has a limiting distibution which is the supremum of the
##' absolute value of a Brownian bridge:
##' \deqn{ p= Pr ( \sup Q \geq q ) = 2 \sum_{i=1}^{\infty} (-1)^{i+1} \exp (-2 i^2 q^2) }{
##'        P(Q >= q) = 2 SUM [i to Inf] (-1)^(i+1).e^(-2.i^2.q^2) }
##' @examples
##' data(kidtran, package="KMsurv")
##' k1 <- kidtran
##' k2 <- k1[k1$gender==1 & k1$race==2, ]
##' c1 <- coxph(Surv(time, delta) ~ age, data = k2)
##' cutp(c1, var="age", plot=TRUE)
##' k2 <- k1[k1$gender==2 & k1$race==2, ]
##' c1 <- coxph(Surv(time, delta) ~ age, data = k2)
##' cutp(c1, var="age")
##' @references Examples are from
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 8.3, pp 273-274.
##' @references Contal C, O'Quigley J, 1999
##' An application of changepoint methods in studying the
##' effect of age on survival in breast cancer.
##' \emph{Computational Statistics & Data Analysis} \bold{30}(3):253--70.
##' \href{http://www.sciencedirect.com/science/article/pii/S0167947398000966}{ScienceDirect}
##'
cutp.coxph <- function(x, ...,
                       var="",
                       plot=FALSE){
    stopifnot(inherits(x, "coxph"))
    stopifnot(!identical(var, ""))
    ## get location to evaluate variables
    ## (i.e. environment if no data frame specified)
    v1 <- get(var, model.frame(x))
### check >2 levels
    stopifnot(length(unique(v1)) > 2)
### get original Surv object
    s1 <- Surv(unclass(model.frame(x)[, 1]))
### need [, 1:2] to remove extra column of all 1's which Surv adds
    s1 <- Surv(s1[,1], s1[,2])
### time, no. events, no. at risk, where deaths occur
    t1 <- tne(s1, onlyEvents=TRUE)
    class(t1) <- "data.frame"
### now check subsets (one for each unique value of variable)
### R1 hold results
    R1 <- cbind(sort(unique(v1)),  NA)
### get log-rank statistic U for each variable
    R1 [, 2] <- sapply(1:nrow(R1), function (j) {
        s2 <- Surv(s1[, "time"][v1 >= R1[j, 1]], s1[, "status"][v1 >= R1[j, 1]])
        t2 <- rbind(c(0, 0, 0), as.data.frame(tne(s2)))
###  fill in first row
        t2[1, ] <- c(0, max(t2$n), 0)
### hold results
        r1 <- rep(0, length=nrow(t1))
### loop is easier to read than apply here
### want to get n and e for each value of time in t1
        for (i in 1:nrow(t1)){
### if time matches, get n
            n1 <- if ( t1$t[i] %in% t2$t ) {
                n1 <- t2$n[which(t2$t == t1$t[i])]
### if time after max t in t2 then is zero
            } else if( t1$t[i] > max(t2$t) ) {
                n1 <- 0
### otherwise get n from:
### n at last preceding time - e at last preceding time
            } else {
                n1 <- tail(t2$n[t2$t < t1$t[i]], 1) - tail(t2$e[t2$t < t1$t[i]], 1)
            }
### if time matches, match no. events from t2,
### otherwise is zero
            e1 <- ifelse(t1$t[i] %in% t2$t, t2$e[which(t2$t == t1$t[i])], 0)
            r1[i] <- e1 - n1*(t1$e[i]/t1$n[i])
        }
### return value
        abs(sum(r1))
    }
                       )
### get cutpoint
    cut1 <- R1[which(R1[, 2] == max(R1[, 2])), ]
###
    findSigma <- function(D){
        (1/(D-1)) * sum(sapply(1:D,
                               ## i in 1:D
                               function(i)
                               ## j in 1:i
                               (1 - sum( sapply( 1:i,
                                                function(j)
                                                1/(D+1-j) ) )) ^2
                               )
                        )
    }
    s2 <- findSigma(sum(t1$e))
    Q1 <- cut1[2] / ( sqrt(s2)*sqrt(sum(t1$e)-1) )
###
    findP <- function(q, acc){
        ## acc = accuracy; should be to Inf but generally 1e3 is enough
        2 * sum(sapply(1:acc, function (j)
                       (-1)^(j+1)* exp(-2*j^2*q^2)
                       )
                )
    }
    p1 <- findP(Q1, 1e3)
    res1 <- data.frame(cp=cut1[1], Q=Q1, p=p1)
    if (plot){
        m1 <- paste0("Test statistic for cut points \n For variable ", var,
                     "\nLarger values indicate cut point more likely here")
        graphics::plot(R1[, 1], R1[, 2],
                       xlab="Cut point",
                       ylab="Test statistic",
                       main=m1,
                       cex=3,
                       ...)
         lines(R1[, 1], R1[, 2])
    }
    return(res1)
}
##' @rdname cutp
##' @aliases cutp.survfit
##' @method cutp survfit
##' @export
##' 
cutp.survfit <- function(x, ..., var="", plot=FALSE){
    f1 <- deparse(x$call)
    f1 <- sub("survfit", "coxph", f1)
    c1 <- eval(parse(text=f1))
    cutp(c1, var=var, plot=plot)
}
