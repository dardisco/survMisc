##' @name tne
##' @title Time, No. at risk, No. events
##' @description For a \code{Surv} object, gives
##' time, no. at risk and no. events.
##' \cr \cr
##' For a \code{survfit} or \code{coxph} object describing right censored data,
##' gives time, no. at risk and no. events by strata.
##' \cr
##' Also (optionally) per-stratum estimates of no. at risk,
##' no. expected and no. events - no. expected.
##' \cr \cr
##' No. events expected (per predictor) is given by:
##' \deqn{\frac{e_i(n[p]_i)}{n_i}}{
##'  e(i)n[p](i) / n(i)}
##' where \eqn{n[p]_i} is the no. at risk for the predictor.
##' @param x A object of class \code{Surv} or \code{survfit}
##' @param ... Additional arguments, e.g. for \code{survfit} object:
##' \cr \cr
##' \code{pred} : if \code{TRUE}
##' gives predicted values as below
##' \cr \cr
##' \code{asList} : if \code{TRUE}
##' will return a list with one element for each stratum
##' @param onlyEvents if \code{TRUE}
##' shows only times at which at least one event occurred.
##' Otherwise shows \emph{all} times recorded (including those censored).
##' @return For a \code{Surv} object: A data frame with columns:
##' \item{t}{time}
##' \item{n}{no. at risk}
##' \item{e}{no. events}
##' \cr \cr
##' For a \code{survfit} or \code{coxph} object:
##' A data frame with columns:
##'  \item{t}{time}
##'  \item{n}{no. at risk (by strata)}
##'  \item{e}{no. events (by strata)}
##'  \item{s}{strata}
##' \cr \cr
##' If \code{asList} = \code{TRUE} then instead a \code{list}
##' with one element for each stratum, where each
##' elements is a \code{data.frame} with columns
##' \strong{t}, \strong{n} and \strong{e} as above.
##' \cr \cr
##' If \code{pred} = \code{TRUE} then instead a \code{data.frame}
##' with columns:
##'  \item{t}{time}
##'  \item{n}{no. at risk (total)}
##'  \item{e}{no. events (total)}
##'  \item{s}{strata}
##'  \item{ns}{no. at risk (by strata)}
##'  \item{Es}{no. events expected (by strata)}
##'  \item{e_Es}{no. events minus no. events expected}
##' @references Example using \code{kidney} data is from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 7.2, pg 210.
##' @rdname tne
##' @export tne
##'
tne <- function(x, ...){
    UseMethod("tne")
    }
##' @rdname tne
##' @aliases tne.Surv
##' @method tne Surv
##' @S3method tne Surv
##' @examples
##' ### Surv object
##' df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
##'                   e=rep(c(0,1),4))
##' s1 <- Surv(df0$t,df0$e,type="right")
##' tne(s1)
##' tne(s1, onlyEvents=FALSE)
##'
tne.Surv <- function(x, ..., onlyEvents=FALSE){
    if(!class(x)=="Surv") stop
    "Only applies to class 'Surv'"
    if(!attr(x, which="type")=="right") warning
    "Only applies to right censored data"
### no. subjects leaving study in each time period
    r1 <- as.matrix(stats::aggregate(x[,"time"],
                                     by=list(x[,"time"]),
                                     FUN=length))
### no. events in each time period
    r2 <- as.matrix(stats::aggregate(x[,"status"],
                                     by=list(x[,"time"]),
                                     FUN=sum))
### join by time
    m1 <- merge(r1,r2,by=1)
    nr1 <- nrow(x)
### make no. at risk
    m1$n <- if(nrow(m1)==1){
        1} else {
            c(nrow(x),nrow(x)-(cumsum(m1[1:(nrow(m1)-1),2])))
        }
### reorder
    m1 <- m1[,c(1,4,3)]
    colnames(m1) <- c("t","n","e")
### (optional) show only times where events occurred
### get ellipsis
    el1 <- list(...)
    if ("onlyEvents" %in% names(el1)){
        onlyEvents <- el1$onlyEvents
    }
    if (onlyEvents) m1 <- m1[!m1[ ,3]==0, ]
    return(m1)
}
###
###----------------------------------------
###
##' @rdname tne
##' @aliases tne.survfit
##' @method tne survfit
##' @S3method tne survfit
##' @examples
##' ### survfit object
##' data(kidney, package="KMsurv")
##' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' tne(s1)
##' tne(s1, asList=TRUE)
##' tne(s1, pred=TRUE)
##' s2 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
##' tne(s2)
##' data(larynx, package="KMsurv")
##' s1 <- survfit(Surv(time, delta) ~ factor(stage) + age, data=larynx)
##' tne(s1, pred=TRUE)
##' data(bmt, package="KMsurv")
##' s1 <- survfit(Surv(t2, d3) ~ z3 +z10, data=bmt)
##' tne(s1)
##' s1 <- survfit(Surv(t2, d3) ~ 1, data=bmt)
##' tne(s1, pred=TRUE, onlyEvents=TRUE)
##'
tne.survfit <- function(x, ..., onlyEvents=FALSE){
    if (!class(x)=="survfit") stop("Only applies to object of class 'survfit'")
###
### create vector for strata identification
    st1 <- NULL
    for(i in 1:length(x$strata)){
### add vector for one strata according to number of rows of strata
        st1 <- c(st1, rep(names(x$strata)[i], x$strata[i]))
    }
### if only one strata (intercept only model)
    if (is.null(x$strata)) st1 <- as.factor(rep(1, length(x$time)))
### create data.table with data from survfit (create column for strata)
    df1 <- data.frame(t=x$time,
                      n=x$n.risk,
                      e=x$n.event,
                      s=factor(st1))
###
### default values for inputs
    pred <- asList <- FALSE
### get ellipsis
    el1 <- list(...)
    if ("pred" %in% names(el1)){
        pred <- el1$pred
    }

    if ("asList" %in% names(el1)){
        asList <- el1$asList
    }
    if (!pred & !asList){
        if (onlyEvents) df1 <- df1[!df1$e==0, ]
        return(df1)
    }
###
###----------------------------------------
###
### get location to evaluate variables
### (i.e. environment or data frame)
    if (is.null(x$call$data)){
        loc1 <- environment(eval(parse(text=as.character(x$call[2]))))
    } else {
        loc1 <- eval(x$call$data)
    }
### length of original data frame (no. rows)
    l1 <- length(get(ls(loc1),loc1))
### hold results:
### time,  event, no. at risk, no. events (total)
### stratum (predictor), no. at risk (per stratum),
### no. Expected events (per stratum)
### no. events - no. Expected events (per stratum)
    df2 <- data.frame(matrix(0, nrow=l1, ncol=7))
    colnames(df2) <- c("t", "n", "e",
                       "s", "ns", "Es", "e_Es")
### get names of time, event and predictor...
### for name of time:
### get formula: eval(x$call[[2]]),
### then get LHS (Surv object): eval(x$call[[2]])[2]
### then remove trailing (): eval(x$call[[2]])[2][[1]]
### then get time variable
    t1 <- as.character(eval(x$call[[2]])[2][[1]][2])
    df2$t <- get(t1, loc1)
    e1 <- as.character(eval(x$call[[2]])[2][[1]][3])
    df2$e <- with(loc1, get(e1))
### strata
     if (!is.null(x$strata)) {
         s1 <- names(x$strata)
         for (i in 1:length(s1)){
### build logical extression to evaluate
### (to check membership of each strata for each row in df2)
             log1 <- unlist(strsplit(s1[i], ","))
             log1 <- sub(" ", "", log1)
             log1 <- unlist(strsplit(log1, "="))
### if formula specified uses factors, then strip them
### e.g. "factor(stage)" --> "stage"
             log1 <- sub("factor", "", log1)
             log1 <- sub("\\(", "", log1)
             log1 <- sub("\\)", "", log1)
             np1 <- length(log1)/2
### ev1 = to evaluate
             ev1 <- NULL
             for (j in seq(1, (np1+1), by=2)){
### expression to evaluate
                 ev2 <- paste("loc1['", log1[j], "'] == ", log1[j+1], sep="")
                 ev1 <- c(ev1, ev2)
             }
             ev1 <- paste0(ev1, collapse=" & ")
### generate expression
             e1 <- parse(text=ev1)
             df2[which(eval(e1)), "s"] <- s1[i]
         }
     } else {
### intercept-only model
         df2$s <- as.factor(rep(1, l1))
     }
    df2 <- df2[order(df2$t, decreasing=FALSE), ]
    rownames(df2) <- NULL
### make no. at risk (total)
    df2[1, "n"] <- nrow(df2)
    df2[2:nrow(df2), "n"] <- nrow(df2)-cumsum(df2[, "e"][-nrow(df2)])
### make no. at risk (per predictor)
    if (!is.null(x$strata)) {
        for (i in 1:length(s1)) {
### total no. at risk
            n1 <- sum(df2$s==s1[i])
            df2[which(df2$s==s1[i]), "ns"] <- seq(n1,1)
        }
    } else {
### intercept-only model
        df2$ns <- df2$n
    }
### make no. expected events (per predictor)
    df2[, "Es"] <- ( df2[, "e"] * df2[, "ns"] ) / df2[, "n"]
### make events - expected
    df2[, "e_Es"] <- df2[, "e"] - df2[, "Es"]
###
    if (!asList){
        if (onlyEvents) df2 <- df2[!df2$e==0, ]
        return(df2)
    }
### list by group
###
### subset and aggregate
    subAgg <- function(i){
        sub1 <- data.frame(t=df2$t[df2$s==s1[i]],
                           ns=df2$ns[df2$s==s1[i]],
                           e=df2$e[df2$s==s1[i]])
        res1 <- as.matrix(
            stats::aggregate(
                sub1, by=list(sub1$t),
                FUN=identity)[2:4])
        if (is.list(res1)){
### no names yet, so refer to columns by no.
### take max for t and np
            res1[,1:2] <- lapply(res1[,1:2], FUN=max)
### take sum for e
            res1[,3] <- lapply(res1[,3], FUN=sum)
            res1 <- matrix(unlist(res1),nrow=nrow(res1),ncol=3)
        }
        colnames(res1) <- c("t", "n", "e")
        if (onlyEvents) res1 <- res1[!res1[, 3]==0, ]
        return(res1)
    }
    res <- sapply(X=1:length(s1), FUN=subAgg)
   ###  names(res) <- paste("pred_", p1, sep="")
    names(res) <- s1
    return(res)
}
##' @rdname tne
##' @aliases tne.coxph
##' @method tne coxph
##' @S3method tne coxph
##' @examples
##' ### coxph object
##' data(kidney, package="KMsurv")
##' c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney)
##' tne(c1)
##' tne(c1, asList=TRUE)
tne.coxph <- function(x, ..., onlyEvents=FALSE){
    f1 <- deparse(x$call)
    f1 <- sub("coxph", "survfit", f1)
    s1 <- eval(parse(text=f1))
###
### default values for inputs
    pred <- asList <- FALSE
### get ellipsis
    el1 <- list(...)
    if ("pred" %in% names(el1)){
        pred <- el1$pred
    }
    if ("asList" %in% names(el1)){
        asList <- el1$asList
    }
    tne(s1, onlyEvents=onlyEvents, pred=pred, asList=asList)
}
