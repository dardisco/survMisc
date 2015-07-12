#' @name tne
#' @title \bold{T}ime, \bold{n}umer at risk and \bold{n}umber of events
#' 
#' @param x A object of class \code{Surv}, \code{survfit},
#' \code{coxph} or \code{formula}.
#' @param ... Additional arguments (not implemented).
#' @param byWhat See \bold{Value} below
#' @param eventsOnly If \code{eventsOnly=TRUE}
#' shows only times at which at least one event occurred.
#' Otherwise shows \emph{all} times recorded
#' (i.e. including those censored)
#' @param shortNames Applies only if \code{byWhat=="list"}
#'  or \code{byWhat=="all"}. The default is to name
#'  the elements of the \code{list} after each stratum.
#'  \cr
#'  As the names for each stratum are made by concatenating the predictor names, this can
#'  become unwieldly.
#'  \cr
#'  If \code{shortNames="FALSE"} they are instead numbered.
#'  \cr
#'  The long names are given as an \code{attribute}
#'  of the \code{data.table}.
#' 
#' @return For a \code{Surv} object: A \code{data.table} with columns:
#'  \item{t}{time}
#'  \item{n}{no. at risk}
#'  \item{e}{no. events}
#' For a \code{survfit}, \code{coxph} or \code{formula}:
#' If \code{byWhat=="time"} (the default), a 
#' \code{data.table} with columns as above. In addition:
#'  \item{s}{stratum; predictor names are separated with an underscore '_'}
#'  \item{ns}{no. at risk (by strata)}
#'  \item{Es}{no. events expected (by strata)}
#'  \item{e_Es}{no. events minus no. events expected}
#' Additional columns returned match those of the predictors in the \code{model.frame}
#' (for \code{survfit} objects) or \code{model.matrix} (in other cases).
#' If \code{byWhat="list"} = then instead a \code{list}
#' with one element for each stratum, where each
#' elements is a \code{data.table} with columns
#' \bold{t}, \bold{n} and \bold{e} as for a \code{Surv} object.
#' If \code{byWhat="all"}, a \code{data.table} with a columns \bold{t}, \bold{n}
#' and \bold{e} as above.
#' There are additional columns for \bold{n} and \bold{e} for each stratum.
#' 
#' @note
#' The number of events expected (per stratum) is given by:
#' \deqn{E = \frac{e_i(n[s]_i)}{n_i}}{
#'       E = e(i)n[s](i) / n(i)}
#' where \eqn{n[s]_i} is the no. at risk for the stratum.
#' \cr \cr
#' If the formula is 'intercept-only', the stratum \code{I=1} is returned.
#' \cr \cr
#' Interaction terms are not currently supported by \code{survfit} objects.
#' @references Example using \code{kidney} data is from:
#' \bold{K&M}. Example 7.2, pg 210.
#' 
#' @rdname tne
#' @export tne
#'
tne <- function(x, ...){
    ## all are methods ultimately passed to tne.data.frame (below)
    ## except tne.Surv()
    UseMethod("tne")
}
#' @rdname tne
#' @aliases tne.Surv
#' @method tne Surv
#' @export
#' 
#' @examples
#' ## Surv object
#' df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
#'                   e=rep(c(0,1),4))
#' s1 <- with(df0, Surv(t, e, type="right"))
#' tne(s1)
#' tne(s1, eventsOnly=FALSE)
#'
tne.Surv <- function(x, ..., eventsOnly=TRUE){
    ## for R CMD check
    n <- e <- status <- NULL
    stopifnot(class(x)=="Surv")
    stopifnot(attributes(x)$type=="right")
    dt1 <- data.table::data.table(unclass(x))
    dt1 <- dt1[, list(n=length(status), e=sum(status)), by=sort(time)]
    dt1 <- dt1[, "n" := c(sum(n), sum(n) - cumsum(n)[ - length(n)])]
    if(eventsOnly) dt1 <- dt1[e==1, ]
    data.table::setnames(dt1, c("t", "n", "e"))
    return(dt1)
}
###
###----------------------------------------
###
#' @rdname tne
#' @aliases tne.coxph
#' @method tne coxph
#' @export
#' 
#' @examples
#' ## coxph object
#' ## K&M 2nd ed. Section 1.2. Table 1.1, page 2.
#' data(hodg, package="KMsurv")
#' hodg <- data.table::data.table(hodg)
data.table::setnames(hodg,
                     c(names(hodg)[!names(hodg) %in% c("score", "wtime")],
                       "Z1", "Z2"))
c1 <- coxph(Surv(time=time, event=delta) ~ Z1 + Z2,
            data=hodg[gtype==1 && dtype==1, ])
tne(c1, byWhat="time")
#' 
#' ## T&G. Section 3.2, pg 47.
#' c1 <- coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc)
#' c1 <- coxph(Surv(time, status==2) ~ edema, data=pbc)
#' c1 <- coxph(Surv(time, status==2) ~ log(bili) + age + cluster(edema), data=pbc)
#' tne(Surv(time, status==2) ~ log(bili) + age, data=pbc)
#' c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney)
#' data(kidney, package="KMsurv")
#' tne(c1)
#' tne(c1, byWhat="row")
#' data(bmt, package="KMsurv")
#' tne(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
#' 
tne.coxph <- function(x, ..., 
                      byWhat=c("time", "row", "cg"),
                      eventsOnly=TRUE,
                      shortNames=FALSE){
    stopifnot(inherits(x, "coxph"))
    byWhat <- match.arg(byWhat)
    x$call$formula <- stats::terms(x=formula(x$call),
                                   specials=c("strata", "cluster", "tt"))
    x$call[[1]] <- as.name("model.frame")
    ## model.frame
    xMF1 <- eval(x$call, parent.frame())
    tne(x=xMF1, ...,
        byWhat=byWhat, 
        eventsOnly=eventsOnly,
        shortNames=shortNames)
}
#' @rdname tne
#' @aliases tne.survfit
#' @method tne survfit
#' @export
#' 
#' @examples
#' ### survfit object
#' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
#' tne(s1)
#' tne(s1, byWhat="stratum")
#' tne(s1, byWhat="time", eventsOnly=TRUE)
#' tne(survfit(Surv(time=time, event=delta) ~ 1, data=kidney))
#' data(larynx, package="KMsurv")
#' tne(survfit(Surv(time, delta) ~ factor(stage) + age, data=larynx))
#' tne(survfit(Surv(t2, d3) ~ z3 +z10, data=bmt), byWhat="all")
#' tne(survfit(Surv(t2, d3) ~ 1, data=bmt))
tne.survfit <- function(x, ...,
                        byWhat=c("time", "row", "cg"),
                        eventsOnly=TRUE,
                        shortNames=FALSE){
    stopifnot(inherits(x, "survfit"))
    byWhat <- match.arg(byWhat)
    c1 <- x$call
    c1[[1]] <- as.name("model.frame")
    mf1 <- eval(c1, parent.frame())
    .getTne(x=mf1, ...,
            byWhat=byWhat, 
            eventsOnly=eventsOnly,
            shortNames=shortNames)
}
###
###----------------------------------------
### 
###
###----------------------------------------
### 
#' @rdname tne
#' @aliases tne.formula
#' @method tne formula
#' @export
#' 
#' @examples
#' ## formula object
#' ## K&M 2nd ed. Example 7.9, pg 224.
#' data(kidney, package="KMsurv")
#' with(kidney, tne(Surv(time=time, event=delta) ~ type, shortNames=TRUE))
#' ## this doesn't work
#' ## s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)
#' tne(Surv(time=t2, event=d3) ~ z3*z10, data=bmt, byWhat="time")
#' tne(Surv(time=t2, event=d3) ~ ., data=bmt)
#' ## example where each list element has only one row
#' ## also names are impractical
#' tne(Surv(time=t2, event=d3) ~ ., data=bmt, byWhat="stratum", shortNames=TRUE)
tne.formula <- function(x, ...,
                      byWhat=c("time", "row", "cg"),
                      eventsOnly=TRUE,
                      shortNames=FALSE){
    byWhat <- match.arg(byWhat)
    stopifnot(inherits(x, "formula"))
### based on code from stats::lm()
    mc1 <- match.call()
    names(mc1)[names(mc1)=="x"] <- "formula"
    mc1 <- mc1[c(1L, match(c("formula", "row"), names(mc1), 0L))]
    mc1$drop.unused.levels <- TRUE
    mc1[[1L]] <- as.name("model.frame")
    mf1 <- eval(mc1, parent.frame())
    tne(x=mf1, ...,
        byWhat=byWhat,
        eventsOnly=eventsOnly,
        shortNames=shortNames)
}
###
###----------------------------------------
###----------------------------------------
### .getTne
### used by tne.coxph, tne.formula and tne.survfit
###----------------------------------------
###----------------------------------------
###
tne.data.frame <- function(x, ...,
                      byWhat=c("time", "row", "cg"),
                      eventsOnly=TRUE,
                      shortNames=FALSE){
    stopifnot(inherits(x, "data.frame"))
    stopifnot(survival::is.Surv(x[[1]]))
    stopifnot(attr(x[[1]], "type") == "right")
### helper functions
    ## paste variables
    ## a function of rows and the data.table
    pasteVar <- function(row1, dt1){
        stopifnot(data.table::is.data.table(dt1))
        paste(
            colnames(dt1)[!colnames(dt1) %in% c("time", "status")],
            row1[!colnames(dt1) %in% c("time", "status")],
            sep="=", collapse=", ")
    }
    ## drop variables
    dropVar <- function(dt1, keep){
        stopifnot(data.table::is.data.table(dt1))
        ## data table columns
        dtCol1 <- colnames(dt1)[!colnames(dt1) %in% keep]
        for(i in seq_along(dtCol1)){
            n1 <- dtCol1[i]
            dt1[, (n1) := NULL]
        }
        return(dt1)
    }
### data.table from x
    xDT <- data.table::data.table(
        cbind(
            stats::model.matrix(terms(x), x),
            stats::model.response(x)))
    xDT <- xDT[order(time)]
#    if(!all(sapply(attr(attr(x, "terms"), "specials"), is.null)))
    ## names of strata/clusters in x
#    xN1 <- grep("^strata\\(.*\\)|cluster\\(.*\\)", names(x), value=TRUE)
    ## names of strata
    xNS1 <- grepl("^strata\\(.*\\)", names(x))
    if(any(xNS1)){
        ## store strata in separate table
        xStr <- data.table::data.table(
            cbind(
                xDT[, list(time, status)],
                x[xNS1]))
        ## remove strata
        xDT <- xDT[, names(xDT)[grepl("^strata\\(.*\\)", names(xDT))] := NULL]
    }
    getTne <- function(dt1, name="str", byWhat=byWhat){
        ## get covariate groups
        dt1[, (name) := as.factor(apply(X=dt1, MARGIN=1, FUN=pasteVar, dt1=dt1))]
        dt1 <- dropVar(dt1=dt1, keep=c("time", "status", name))
        ## number at risk
        dt1[, "n" := rev(seq.int(nrow(x)))]
        ## number at risk per covariate group
        ncg1 <- paste0("n", name)
        dt1[, (ncg1) := rev(seq.int(length(n))), by=eval((name))]
        data.table::setnames(dt1,
                             c("t", "e", name, "n", ncg1))
        data.table::setcolorder(dt1,
                                c("t", "n", "e", name, ncg1))
        ## get long names
        ln1 <- data.table::data.table(
            "id" = dt1[, seq_along(levels(eval(as.name(name))))],
            "name" = dt1[, levels(eval(as.name(name)))])
        if(shortNames){
            dt1[, (name) := as.numeric(eval(as.name(name)))]
            attr(dt1, "longNames") <- ln1
        }
        l1 <- by(dt1[, list(t, n, e, eval(as.name(ncg1)))],
                 dt1[, eval(as.name(name))],
                 identity)
        dimnames(l1) <- list(dt1[, unique(eval(as.name(name)))])
        if(shortNames) attr(l1, "longNames") <- ln1
        if(byWhat="time"){
            dt2 <- data.table::data.table(t=sort(unique(dt1[, t])))
            n1 <- dt1[, unique(eval(as.name(name)))]
            t1 <- dt1[, sort(unique(t))]
            dt2 <- data.table::data.table(
                matrix(rep(0L, 2 * length(n1) * length(t1)),
                       ncol=2 * length(n1),
                       nrow=length(t1)))
#            dt2[, "t" := t1]
            for(k in seq_along(n1)){
                ## list element per covariate group
                el1 <- l1[k]
                ## index (match time)
                ind1 <- as.integer(which(t1 %in% el1[[1]]$t))
                ## number at risk
                if(length(ind1)==1L){
                    set(dt2, i=NULL, j=k, value=el1[[1]]$n)
                } else {
                    set(dt2, i=ind1, j=k, value=el1[[1]]$n)
                for (l in seq_along)
                set(dt2, i=seq.int(ind1[l]), j=k, value=el1[[1]]$n[l])
                
                
                i1 <- dt1[, which(cg==levels(cg)[i])]
            ## n = number at risk
            n1 <- paste0("n_", dt1[, levels(cg)[i]])
            tn1 <- dt1[i1, max(ncg), by=t]
            dt2[t %in% tn1[, t], i := tn1[, V1]]
            dt2[, (n1) := zoo::na.locf(
                           eval(bquote(.(as.name(n1)))),
                           na.rm=FALSE, fromLast=FALSE)]
            dt2[, (n1) := zoo::na.locf(
                           eval(bquote(.(as.name(n1)))),
                           na.rm=FALSE, fromLast=TRUE)]
            ## e = number of events
            e1 <- paste0("e_", dt1[, levels(cg)[i]])
            te1 <- dt1[i1, sum(e, na.rm=TRUE), by=t]
            dt2[t %in% te1[, t], (e1) := te1[, V1]]
            dt2[is.na(eval(bquote(.(as.name(e1))))), (e1) := 0]
        }
        ## make no. at risk (total) per time period
        dt2[, "n" := rowSums(.SD), .SDcols = grep("n_", colnames(dt2))]
        ## total events per time period
        dt2[, "e" := rowSums(.SD), .SDcols = grep("e_", colnames(dt2))]
        data.table::setcolorder(dt2,
                                c("t", "n", "e",
                                  colnames(dt2)[!colnames(dt2) %in% c("t", "n", "e")]))
        return(dt2)
 
        ## 
        l1 <- by(dt1[, list(t, n, e, eval(as.name(ncg1)))],
                 dt1[, eval(as.name(name))],
                 identity)
        dimnames(l1) <- list(dt1[, unique(eval(as.name(name)))])
        if(shortNames) attr(l1, "longNames") <- ln1
        l1
    }

        
        dt[, list(t,n,ncg), ][, list(list(.SD)), by=cg]$V1
        dt[, list(t,n,ncg), by=unique(cg)][, list(list(.SD))]$V1
        [, list(list(.SD)), by=cg]$V1
        dt[, list(t, n, e), mult="first"]
        
        l1 <-           dt[, list(t, "n"=eval(as.name(ncg1)), e, eval(as.name(name)))][, list(list(.SD)), by=eval((name))]$V1
        dt[, list(t, "n"=nstr, e, str)][, list(list(.SD)), by=str]$V1

        l1 <-
            dt[, list(t, "n"=ncg, e, cg) ][, list(list(.SD)), by=cg]$V1

        if(byWhat %in% c("row", "cg") && eventsOnly) dt1 <- dt1[e==1, ]
        
        l1 <- dt[, list(t, "n"=ncg, e, str) ][, list(list(.SD)), by=str]$V1
        
    


    ## get covariate groups
    dt1[, "cg" := as.factor(apply(X=dt1, MARGIN=1, FUN=pasteVar, dt=dt1))]
    dt1 <- dropVar(dt=dt1, keep=c("cg", "time", "status"))
    data.table::setnames(dt1,
                         c("t", "e", "cg", "n"))
    data.table::setcolorder(dt1,
                            c("t", "n", "e", "cg"))
    ## get long names
    ln1 <- data.table::data.table(
        "id" = dt1[, seq_along(levels(cg))],
        "name" = dt1[, levels(cg)])
    if(byWhat %in% c("row", "cg") && eventsOnly) dt1 <- dt1[e==1, ]
###
    if(byWhat == "row"){
        ## expected events
        dt1[, "Ecg" := (e * ncg) / n]
        ## make events - expected
        dt1[, "e_Ecg" := e - Ecg]
        if (shortNames){
            dt1[, cg := as.numeric(cg)]
            attr(dt1, "longNames") <- ln1
        }
        return(dt1)
    }
###
    if(byWhat == "time"){
        if(exists("xStr")){
            dt1[, "str" := xStr[, str]]
            
   

            
        dt1[, "st" := NULL]
        ### slow!
        tneTime <- function(dt){
            dt1 <- dt
        dt2 <- data.table::data.table(t=sort(unique(dt1[, t])))
        for(i in dt1[, 1:length(levels(cg))]){
            ## index (match covariate group)
            i1 <- dt1[, which(cg==levels(cg)[i])]
            ## n = number at risk
            n1 <- paste0("n_", dt1[, levels(cg)[i]])
            tn1 <- dt1[i1, max(ncg), by=t]
            dt2[t %in% tn1[, t], (n1) := tn1[, V1]]
            dt2[, (n1) := zoo::na.locf(
                           eval(bquote(.(as.name(n1)))),
                           na.rm=FALSE, fromLast=FALSE)]
            dt2[, (n1) := zoo::na.locf(
                           eval(bquote(.(as.name(n1)))),
                           na.rm=FALSE, fromLast=TRUE)]
            ## e = number of events
            e1 <- paste0("e_", dt1[, levels(cg)[i]])
            te1 <- dt1[i1, sum(e, na.rm=TRUE), by=t]
            dt2[t %in% te1[, t], (e1) := te1[, V1]]
            dt2[is.na(eval(bquote(.(as.name(e1))))), (e1) := 0]
        }
        ## make no. at risk (total) per time period
        dt2[, "n" := rowSums(.SD), .SDcols = grep("n_", colnames(dt2))]
        ## total events per time period
        dt2[, "e" := rowSums(.SD), .SDcols = grep("e_", colnames(dt2))]
        data.table::setcolorder(dt2,
                                c("t", "n", "e",
                                  colnames(dt2)[!colnames(dt2) %in% c("t", "n", "e")]))
        return(dt2)
        }
        if(eventsOnly) dt2 <- dt2[e >= 1, ]
        if (shortNames){
            sn1 <- dt1[, seq_along(levels(cg))]
            ne1 <- c("n_", "e_")
            n1 <- as.vector(outer(ne1, sn1, paste, sep=""))
            data.table::setnames(dt2, c("t", "n", "e", n1))
            attr(dt2, "longNames") <- ln1
        }
        return(dt2)
    }
### 
    if(byWhat == "cg"){
        ## make a list
        l1 <- dt1[, list(t, "n"=ns, e, cg) ][, list(list(.SD)), by=cg]$V1
        if(shortNames){
            attr(l1, "longNames") <- ln1
        } else {
            names(l1) <- dt1[, levels(cg)]
        }
        return(l1)
    }
}




 
