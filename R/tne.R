##' @name tne
##' @title Time, No. at risk, No. events
##' 
##' @param x A object of class \code{Surv}, \code{survfit},
##' \code{coxph} or \code{formula}.
##' @param ... Additional arguments (not implemented)
##' @param eventsOnly If \code{eventsOnly=TRUE}
##' shows only times at which at least one event occurred.
##' Otherwise shows \emph{all} times recorded
##' (i.e. including those censored)
##' @param what See \bold{Value} below
##' @param nameStrata Applies only if \code{what=="list"}
##' or \code{what=="all"}. The default is to name
##' the elements of the \code{list} after each stratum.
##' \cr
##' As the names for each stratum are made by concatenating the predictor names, this can
##' become unwieldly.
##' \cr
##' If \code{nameStrata="FALSE"} they are instead numbered. A list is returned
##' with the numbered \code{list} or \code{data.table} 
##' and a \code{vector} giving the names of the strata.
##' 
##' @return For a \code{Surv} object: A \code{data.table} with columns:
##' \item{t}{time}
##' \item{n}{no. at risk}
##' \item{e}{no. events}
##' For a \code{survfit}, \code{coxph} or \code{formula}:
##' If \code{what="table"} (the default), a 
##' \code{data.table} with columns as above. In addition:
##'  \item{s}{stratum; predictor names are separated with an underscore '_'}
##'  \item{ns}{no. at risk (by strata)}
##'  \item{Es}{no. events expected (by strata)}
##'  \item{e_Es}{no. events minus no. events expected}
##' Additional columns returned match those of the predictors in the \code{model.frame}
##' (for \code{survfit} objects) or \code{model.matrix} (in other cases).
##' If \code{what="list"} = then instead a \code{list}
##' with one element for each stratum, where each
##' elements is a \code{data.table} with columns
##' \bold{t}, \bold{n} and \bold{e} as for a \code{Surv} object.
##' If \code{what="all"}, a \code{data.table} with a columns \bold{t}, \bold{n}
##' and \bold{e} as above.
##' There are additional columns for \bold{n} and \bold{e} for each stratum.
##' 
##' @note
##' The number of events expected (per stratum) is given by:
##' \deqn{E = \frac{e_i(n[s]_i)}{n_i}}{
##'       E = e(i)n[s](i) / n(i)}
##' where \eqn{n[s]_i} is the no. at risk for the stratum.
##' \cr \cr
##' If the formula is 'intercept-only', the stratum \code{I=1} is returned.
##' \cr \cr
##' Interaction terms are not currently supported by \code{survfit} objects.
##' @references Example using \code{kidney} data is from:
##' \bold{K&M}. Example 7.2, pg 210.
##' 
##' @rdname tne
##' @export tne
##'
tne <- function(x, ...){
    UseMethod("tne")
### all are methods ultimately passed to .getTne (below)
### except tne.Surv()
}
##' @rdname tne
##' @aliases tne.Surv
##' @method tne Surv
##' @export
##' 
##' @examples
##' ### Surv object
##' df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
##'                   e=rep(c(0,1),4))
##' s1 <- Surv(df0$t, df0$e, type="right")
##' tne(s1)
##' tne(s1, eventsOnly=TRUE)
##'
tne.Surv <- function(x, ..., eventsOnly=FALSE){
### for R CHD check
    n <- e <- status <- NULL
    if(!class(x)=="Surv") stop
    "Only applies to class 'Surv'"
    if(!attr(x, which="type")=="right") warning
    "Only applies to right censored data"
    dt1 <- data.table(unclass(x))
    dt1 <- dt1[, list(n=length(status), e=sum(status)), by=sort(time)]
    dt1 <- dt1[, "n" := c(sum(n), sum(n) - cumsum(n)[ - length(n)])]
    if(eventsOnly) dt1 <- dt1[e==1, ]
    setnames(dt1, c("t", "n", "e"))
    return(dt1)
}
###
###----------------------------------------
###
##' @rdname tne
##' @aliases tne.survfit
##' @method tne survfit
##' @export
##' 
##' @examples
##' ### survfit object
##' data(kidney, package="KMsurv")
##' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' tne(s1)
##' tne(s1, what="all")
##' tne(s1, what="all", eventsOnly=TRUE)
##' tne(survfit(Surv(time=time, event=delta) ~ 1, data=kidney))
##' data(larynx, package="KMsurv")
##' tne(survfit(Surv(time, delta) ~ factor(stage) + age, data=larynx))
##' data(bmt, package="KMsurv")
##' tne(survfit(Surv(t2, d3) ~ z3 +z10, data=bmt), what="all")
##' tne(survfit(Surv(t2, d3) ~ 1, data=bmt))
tne.survfit <- function(x, ...,
                        eventsOnly=FALSE,
                        what=c("table", "list", "all"),
                        nameStrata=TRUE){
    stopifnot(class(x)=="survfit")
    what <- match.arg(what) 
    mf <- x$call
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
### model terms
    mt <- attr(mf, "terms")
    stopifnot(is.empty.model(mt)==FALSE)
### make data.table from model frame excluding response
    dt1 <- data.table(mf[-1])
### check if intercept-only model
    s1 <- attr(mt, "term.labels")
    if (length(s1)==0) {
        dt1 <- data.table("I" = rep(1, nrow(mf)))
    }
    return( .getTne(dt1, mf=mf, eventsOnly=eventsOnly, what=what) )
}
###
###----------------------------------------
### 
##' @rdname tne
##' @aliases tne.coxph
##' @method tne coxph
##' @export
##' 
##' @examples
##' ### coxph object
##' data(kidney, package="KMsurv")
##' c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney)
##' tne(c1)
##' tne(c1, what="list")
##' tne(coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
##' 
tne.coxph <- function(x, ...,
                      eventsOnly=FALSE,
                      what=c("table", "list", "all"),
                      nameStrata=TRUE){
    stopifnot(class(x)=="coxph")
    what <- match.arg(what) 
    mf <- x$call
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
### model terms
    mt <- attr(mf, "terms")
    stopifnot(is.empty.model(mt)==FALSE)
### if some terms, exclude intercept
    if (length( attr(mt, "term.labels") ) > 0) {
        attr(mt, "intercept") <- 0
    }
### get model matrix
    dt1 <- data.table(model.matrix(mt, mf, contrasts))
    return(.getTne(dt1, mf, eventsOnly=eventsOnly, what=what, nameStrata=nameStrata))
}
###
###----------------------------------------
### 
##' @rdname tne
##' @aliases tne.formula
##' @method tne formula
##' @export
##' 
##' @examples
##' ### formula object
##' data(kidney, package="KMsurv")
##' ### this doesn't work
##' ### s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)
##' tne(Surv(time=t2, event=d3) ~ z3*z10, data=bmt, what="all")
##' tne(Surv(time=t2, event=d3) ~ ., data=bmt)
##' ### example where each list element has only one row
##' ### also names are impractical
##' tne(Surv(time=t2, event=d3) ~ ., data=bmt, what="list", nameStrata=FALSE)
tne.formula <- function(x, ...,
                        eventsOnly=FALSE,
                        what=c("table", "list", "all"),
                        nameStrata=TRUE){
    what <- match.arg(what) 
### code copied from lm()
    mf <- match.call()
    if (class(x)=="formula") names(mf)[names(mf)=="x"] <- "formula"
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
### model terms
    mt <- attr(mf, "terms")
    stopifnot(is.empty.model(mt)==FALSE)
### if some terms, exclude intercept
    if (length( attr(mt, "term.labels") ) > 0) {
        attr(mt, "intercept") <- 0
    }
### get model matrix
    dt1 <- data.table(model.matrix(mt, mf, contrasts))
    return( .getTne(dt1, mf, eventsOnly=eventsOnly,
                    what=what, nameStrata=nameStrata) )
}
###
###----------------------------------------
###----------------------------------------
### .getTne
### used by tne.coxph, tne.formula and tne.survfit
###----------------------------------------
###----------------------------------------
###
.getTne <- function(dt1, mf,
                    eventsOnly=FALSE,
                    what,
                    nameStrata=TRUE){
    stopifnot(what %in% c("table", "list", "all"))
### 
### 'dt1' is a data.table as returned by one of the functions above
### It contains the model.matrix or model.frame (excluding response) for the object
### It excludes the intercept term, if other terms are present
### 
### 'mf' is the model frame for the object
### 
### for R CMD check
    s <- Es <- .SD <- NULL
    n1 <- names(dt1)
### 
### build up expression to evalate in dt1
### get names of predictors
### use qQuote to preserve quotes in quotes
### i.e. change " to \"
### avoid fancy quotes in UNIX!
### (note this is typically recommended only for display to terminal)
### 
    fq1 <- unname(unlist(options("useFancyQuotes")))
    options(useFancyQuotes=FALSE)
    f1 <- function(x) paste(dQuote(paste(x, "=", sep="")),
                                   ", get(", dQuote(x), ")", sep="")

    t1 <- paste(sapply(n1, f1), sep=",")
### if more than one predictor
### add space before each name (except for first name)
    if(length(t1) >1 ) {
        t1[2:length(t1)] <- sub("\"", "\"_", t1[2:length(t1)], fixed=TRUE)
    }
    t1 <- paste(t1, collapse=",")
    t2 <- paste("paste(", t1, ", sep='')" )
    p1 <- parse(text=t2)
    q <- quote(eval(p1))
    options(useFancyQuotes=fq1)
### make one column indicating strata
### (one value for each combination of predictors)
### on Linux, argument envir=.SD not necessary
### include SD (subset data.table) to allow eval to work there
### dt1[, "s" := as.factor(eval(q, envir=.SD))]
    dt1[, "s" := dt1[, eval(q, envir=.SD)]]
    dt1[, "s" := as.factor(s)]
###
    stopifnot(attr(model.response(mf), "type")=="right")
    y <- data.table(unclass(model.response(mf, "numeric")))
    dt1[, c("t", "e") := y]
    dt1 <- dt1[order(t)]
### number at risk
    dt1[, "n" := c(nrow(dt1), nrow(dt1)- cumsum(e)[-nrow(dt1)]) ]
### number at risk per strata
    dt1[, "ns" := rev(seq(length(n))), by=s]
    nc1 <- ncol(dt1)
    setcolorder(dt1,
                c(nc1-3, nc1-1, nc1-2, nc1, nc1-4, 1:(nc1-5)))
###
###----------------------------------------
### table
###----------------------------------------
### 
    if(what=="table"){
### otherwise
### ### make no. expected events (per predictor)
        if(eventsOnly) {
            dt1 <- dt1[e==1, ]
        }
        dt1[, "Es" := (e * ns) / n]
### ### make events - expected
        dt1[, "e_Es" := e - Es ]
        nc1 <- ncol(dt1)
        setcolorder(dt1,
                    c(1:3, 5, 4, nc1-1, nc1, 6:(nc1-2))
                    )
        return(dt1)
    }
### 
###----------------------------------------
### list
###----------------------------------------
### 
### make list, one for each stratum
### need call to .SD to make 2nd list work
    l1 <- dt1[, list(t, "n"=ns, e, s) ][, list(list(.SD)), by=s]$V1
    if(what=="list") {
        if(eventsOnly) {
            dt1 <- dt1[e==1, ]
        }
        if(nameStrata) {
            names(l1) <- unique(dt1$s)
            return(l1)
        } else {
            return(list(strata=dt1$s, tne=l1))
        }
    }
### 
###----------------------------------------
### else what=="all"
### i.e. table with one column per stratum
###----------------------------------------
### 
    for (i in seq_along(l1)){
### ### get no. at risk and no. events per time
        l1[[i]] <- l1[[i]][, list(n=max(n), e=sum(e)), by=t]
### ### rename 'n' and 'e' columns to avoid duplicate column names
### ### needed if merging >3 data.frames
        setnames(l1[[i]], c("t", paste0("n", i), paste0("e", i)))
    }
### merge all elements in list (rather inefficient)
### need allow.cartesian to prevent error
### if new no. rows > no. rows in longest element in list
    m1 <- data.table::data.table(Reduce(function(...)
                                        merge(..., by="t",
                                              all=TRUE, allow.cartesian=TRUE),
                                        l1))
### for 'n' carry last observation back
### (to fill in missing values in first rows)
    data.table::set(m1, j=grep("n", colnames(m1)),
                    value=zoo::na.locf(m1[, grep("n", colnames(m1)),
                    with=FALSE],
                    fromLast=TRUE))
### for remaining 'n' and 'e', replace NA with zero
### for 'n' this will be elements in the tail of the vector
    for (j in seq_len(ncol(m1))){
        data.table::set(m1, which(is.na(m1[[j]])), j, as.integer(0))
    }
### initialise new constants (for R CMD check)
    n <- e <- NULL
### make no. at risk (total) per time period
    m1[, n := rowSums(.SD), .SDcols = grep("n", colnames(m1))]
### total events per time period
    m1[, e := rowSums(.SD), .SDcols = grep("e", colnames(m1))]
    setcolorder(m1,
                c(1, ncol(m1) - 1, ncol(m1), 2:(ncol(m1)-2))
                )
###
    if(eventsOnly) {
        m1 <- m1[e >= 1, ]
    }
### 
###----------------------------------------
### change names if required
###----------------------------------------
### 
    if (nameStrata){
        s1 <- levels(dt1$s)
        n1 <- c("n_", "e_")
        n1 <- as.vector(outer(n1, s1, paste, sep=""))
        data.table::setnames(m1, c("t", "n", "e", n1))
        return(m1)
    }
    return(list(strata=s1, data=m1))
}
