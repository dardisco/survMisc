#' @name tn
#' @title \bold{t}ime and \bold{n}umber at risk, with status or
#' numer of events.
#'
#' @include printTn.R
#'
#' @param x
#' For the default method, a \code{numeric} vector indicating
#' \emph{status}.
#'  \cr
#' Each element indicates whether an event occurred (\code{1}) or
#' not (\code{0}) for an observation.
#'  \cr
#' These are assumed to be ordered by discrete times.
#'  \cr
#' This is similar to the \code{event} argument for \code{Surv}
#' objects.
#'  \cr
#' Methods are available for objects of class
#' \code{Surv}, \code{survfit},
#' \code{coxph} and \code{formula}.
#'  \cr
#' The methods for \code{data.frame} (for a model frame)
#' and \code{data.table} are not intended for interactive use.
#'  \cr
#' The method for \code{tn} objects converts the object
#' from 'long' to 'wide' form or \emph{vice versa}.
#' @param shape See \bold{Value} below.
#' @param tweOnly \bold{T}imes \bold{w}ith \bold{e}vnets.
#'  \cr
#' If \code{tweOnly=TRUE} (the default),
#' show only times at which \emph{at least one} event occurred.
#' Otherwise shows \emph{all} times recorded
#' (i.e. including those censored).
#' @param abbNames
#' \bold{Abb}reviate names.
#'  \cr
#' The default is to name
#' the elements of the \code{list} after each stratum.
#'  \cr
#'  As the names for each stratum are made by concatenating
#'  the predictor names, this can
#'  become unwieldly.
#'  \cr
#'  If \code{abbNames="FALSE"} they are instead numbered.
#'  \cr
#' In either case, the \code{longNames} are given
#' as an \code{attribute} of the returned \code{tn} object.
#' @param contrasts.arg Methods for handling factors.
#'  \cr
#' A \code{list}. The \code{names} are the names of
#' columns of the \code{model.frame} containing
#' \code{factor}s.
#'  \cr
#' The \emph{values} are used as replacement
#' values for the \code{contrasts} replacement function.
#' These should be functions (given as character strings)
#' or numeric matrices. 
#'  \cr
#' This can be passed from
#' \code{survfit}, \code{coxph} and \code{formula} objects to:
#'  \cr
#' ?stats::model.matrix
#' @param call Used to pass the \code{call} from a \code{formula}
#'  to the final \code{tn.data.table} method. 
#'
#' @return A \code{data.table} with the additional \code{class}
#'  \code{tn}. The same information can be displayed in two ways.
#'  \cr
#' If \code{shape=="wide"} (the default), the data is
#'  returned in 'wide' format. There is one
#'  row for each time.
#'  \cr
#'  For a \code{numeric} or \code{Surv} object this has columns: 
#'   \item{t}{time.}
#'   \item{n}{number at risk.}
#'   \item{e}{Number of events.}
#'  A \code{survfit}, \code{coxph} or \code{formula} object
#'  will have additional columns for \code{n} and \code{e}
#'  for each covariate group. 
#'
#' If \code{shape=="obs"}, the data is
#'  returned in 'long' format. There is one
#'  row for each observation of an event.
#'  \cr
#'  For a \code{numeric} or \code{Surv} object this has columns: 
#'   \item{time}{time}
#'   \item{n}{number at risk}
#'   \item{status}{event observed (\code{=1}) or not \code{=0}.}
#'  A \code{survfit}, \code{coxph} or \code{formula} object
#'  will have additional columns:
#'   \item{cg}{\bold{c}ovariate \bold{g}roup.
#'  \cr
#'  This is formed from combining the variables; these
#'  are separated by an comma ','}
#'   \item{ncg}{\bold{n}umber at risk, by \bold{c}ovariate \bold{g}roup}
#'
#' \bold{Special terms}.
#'  \cr
#'  The following are considered 'special'
#'  terms in a survival model:
#'   \item{strata}{For a stratified model, a \code{by list} with
#'    one \code{tn} object per strata as above.}
#'   \item{cluster}{These terms are dropped.}
#'   \item{tt}{The variable is unchhanged. That is, time-transform
#'    terms are handled as if the the function
#'    \code{tt(x)} was \code{identity(x)}.
#' \cr
#' \bold{Attribures}.
#'  \cr
#'  The returned object will also have attributes for
#'  the following:
#'   \item{shape}{as above}
#'   \item{tweOnly}{as above}
#'   \item{abbNames}{as above}
#'   \item{ncg}{Number of covariate groups}
#'   \item{call}{The call used to generate the object}
#'
#' @note Partial matching is available for the
#' following arguments, based on the characters in bold:
#'  \cr
#' \bold{sh}ape, \bold{abb}Names, \bold{twe}Only,
#' \bold{con}trasts.arg
#'  \cr
#' Currently only binary status and right-censoring
#' are supported.
#'
#' @rdname tn
#' @export
#'
## for R CMD check
n <- e <- status <- NULL
tn <- function(x, ...) UseMethod("tn")
## all are methods ultimately passed to
## tn.data.frame (below)
## except tn.numeric() and tn.Surv()
###----------------------------------------
#' @rdname tn
#' @aliases tn.numeric
#' @method tn numeric
#' @export
#'
tn.numeric <- function(x, ...,
                       tweOnly=TRUE){
    partMatch(env1=environemt(), ...)
    stopifnot(all(x %in% c(0,1)))
    res1 <- data.table::data.table(
        "t"=(t <- seq_along(x)),
        "status"=x,
        "n"=rev(t))
    if(tweOnly) res1 <- res1[status==1, ]
    setAttrTn(res1,
              shape="long",
              tweOnly=tweOnly,
              abbNames=TRUE,
              ncg=0,
              call=match.call(),
              class=c("tn", class(res1)))
    return(res1)
}
#' @rdname tn
#' @aliases tn.Surv
#' @method tn Surv
#' @export
#'
tn.Surv <- function(x, ..., 
                    shape=c("wide", "long"),
                    tweOnly=TRUE,
                    call=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environemt(), ...)
    if(is.null(call)) call <- match.call()
    stopifnot(inherits(x, "Surv"))
    stopifnot(attributes(x)$type=="right")
    res1 <- data.table::data.table(unclass(x))
    if(shape=="long"){
        res1[, "n" := rev(seq.int(nrow(res1)))]
        data.table::setcolorder(res1, c("time", "n", "status"))
        data.table::setnames(res1, c("t", "n", "status"))
        if(tweOnly){
            data.table::setkey(res1, status)
            res1 <- res1[status > 0.5, ]
        }
    }
    if(shape=="wide"){
        res1 <- res1[, list("n"=length(status),
                            "e"=sum(status)),
                     by=sort(time, na.last=TRUE)]
        res1[, "n" := c(sum(n), sum(n) - cumsum(n)[ - length(n)])]
        data.table::setnames(res1, c("t", "n", "e"))
        if(tweOnly){
            data.table::setkey(res1, e)
            res1 <- res1[e > 0.5, ]
        }
    }
    ## sort by time or t
    data.table::setkeyv(res1, colnames(res1)[1])
    setAttrTn(res1,
              shape=shape,
              tweOnly=tweOnly,
              ncg=0,
              call=call,
              class=c("tn", class(res1)))
    return(res1)
}
###
###----------------------------------------
###
#' @rdname tn
#' @aliases tn.coxph
#' @method tn coxph
#' @export
#'
tn.coxph <- function(x, ...,
                     shape=c("wide", "long"),
                     tweOnly=TRUE,
                     abbNames=TRUE,
                     contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    x$call$formula <- stats::terms(
        x=formula(x$call),
        specials=c("strata", "cluster", "tt"))
    x$call$drop.unused.levels <- TRUE
    call1 <- x$call
    x$call[[1]] <- as.name("model.frame")
    ## model.frame
    xMF1 <- eval(x$call, parent.frame())
    tn(x=xMF1,
       shape=shape,
       tweOnly=tweOnly,
       abbNames=abbNames,
       contrasts.arg=contrasts.arg,
       call=call1)
}
#' @rdname tn
#' @aliases tn.survfit
#' @method tn survfit
#' @export
#'
tn.survfit <- function(x, ..., 
                       shape=c("wide", "long"),
                       tweOnly=TRUE,
                       abbNames=TRUE,
                       contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    x$call$formula <- stats::terms(
        x=formula(x$call),
        specials=c("strata", "cluster", "tt"))
    x$call$drop.unused.levels <- TRUE
    call1 <- x$call
    x$call[[1]] <- as.name("model.frame")
    xMF1 <- eval(x$call, parent.frame())
    tn(x=xMF1,
       shape=shape,
       tweOnly=tweOnly,
       abbNames=abbNames,
       contrasts.arg=contrasts.arg,
       call=call1)
}
#' @rdname tn
#' @aliases tn.formula
#' @method tn formula
#' @export
#'
tn.formula <- function(x, ...,
                       shape=c("wide", "long"),
                       tweOnly=TRUE,
                       abbNames=TRUE,
                       contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    stopifnot(inherits(x, "formula"))
### based on code from stats::lm()
    mc1 <- match.call()
    names(mc1)[names(mc1)=="x"] <- "formula"
    mc1 <- mc1[c(1L, match(c("formula", "data"), names(mc1), 0L))]
    mc1$drop.unused.levels <- TRUE
    call1 <- mc1
    mc1[[1L]] <- as.name("model.frame")
    mf1 <- eval(mc1, parent.frame())
    tn(x=mf1,
       shape=shape,
       tweOnly=tweOnly,
       abbNames=abbNames,
       contrasts.arg=contrasts.arg,
       call=call1)
}
#' @rdname tn
#' @aliases tn.data.frame
#' @method tn data.frame
#' @export
#'
tn.data.frame <- function(x, ...,
                          shape=c("wide", "long"),
                          tweOnly=TRUE,
                          abbNames=TRUE,
                          contrasts.arg=NULL,
                          call=NULL){
    stopifnot(survival::is.Surv(x[[1]]))
    stopifnot(attr(x[[1]], "type") == "right")
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
### data.table from x
    xDT <- data.table::data.table(
        cbind(
            stats::model.matrix(terms(x), x,
                                contrasts.arg=contrasts.arg),
            stats::model.response(x)))
### data.table::setkey(xDT, time)
    ## names of clusters
    xNC1 <- grepl("^cluster\\(.*\\)", names(x))
    if(any(xNC1)){
        ## drop cluster terms
        xDT[, names(xDT)[grepl("^cluster\\(.*\\)",
                               names(xDT))] := NULL]
    }
    ## names of strata
    xNS1 <- grep("^strata\\(.*\\)", names(x))
    if(any(xNS1)){
        ## strata numbers
        xDTstn1 <- grep("^strata\\(.*\\)", names(xDT))
        ## separate table onlly for strata
        xDTstr <- xDT[, .SD, .SDcols=xDTstn1]
        collapseDT(xDTstr, except=NA, nName="strat")
        data.table::set(xDT, j=xDTstn1, value=NULL)
        collapseDT(xDT,
                   except=c("time", "status"),
                   nName="cg")
        xDT[, "cg" := as.factor(cg)]
        xDT[, "strat" := as.factor(xDTstr[, strat])]
        l1 <- by(xDT[, list(time, status, cg)],
                 xDT[, strat],
                 identity)
        res1 <- vector(mode="list")
        for (i in seq_along(l1)){
            res1[i] <- tn(l1[[i]], 
                          shape=shape,
                          tweOnly=tweOnly,
                          abbNames=abbNames)
            names(res1)[i] <- names(l1)[i]
        }
        data.table::setattr(res1,
                            "class",
                            c("stratTn", class(res1)))
        data.table::setattr(res1, "call", call)
        return(res1)
    } else {
        if(stats::is.empty.model(x)){
            ## convert to Surv object
            s1 <- xDT[, Surv(time, status)]
            return(tn(s1,
                      shape=shape,
                      tweOnly=tweOnly,
                      call=call))
        } else {
            collapseDT(xDT,
                       except=c("time", "status"),
                       nName="cg")
            tn(x=xDT,
               shape=shape,
               tweOnly=tweOnly,
               abbNames=abbNames,
               call=call)
        }
    }
}
#' @rdname tn
#' @aliases tn.data.table
#' @method tn data.table
#' @export
#' 
tn.data.table <- function(x, ...,
                          shape=c("wide", "long"),
                          tweOnly=TRUE,
                          abbNames=TRUE,
                          call=NULL){
    stopifnot(all(names(x) %in% c("time", "status", "cg")))
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    data.table::setkey(x, time, cg)
    ## number at risk
    x[, "n" := rev(seq.int(nrow(x)))]
    ## number at risk per covariate group
    x[, "ncg" := rev(seq.int(length(n))), by=cg]
    ## get long names
    x[, "cg" := as.factor(cg)]
    ln1 <- data.table::data.table(
        "id" = x[, unique(as.integer(cg))],
        "longName" = x[, levels(cg)])
    if(abbNames) x[, "cg" := as.integer(cg)]
    if(tweOnly){
        ## times with greater than zero events
        x[, "tGr0e" := sum(status) > 0L, by=time]
        ## subset longNames
        cg1 <- x[x[, tGr0e], unique(as.integer(cg))]
        ln1[, "twe" := id %in% cg1]
        ## subset x and drop tGr0e from x
        x <- x[(tGr0e)]
        x[, tGr0e := NULL]
    }
    data.table::setnames(x, c("t", colnames(x)[-1]))
    data.table::setcolorder(x,
                            c("t", "n", "status", "cg", "ncg"))
    setAttrTn(x,
              "shape"="long",
              "tweOnly"=tweOnly,
              "abbNames"=abbNames,
              "longNames"=ln1,
              "ncg"=nrow(ln1),
              "call"=call,
              "class"=c("tn", class(x)))
    if(shape=="long") return(x)
    tn(x=x,
       shape="long",
       tweOnly=tweOnly,
       abbNames=abbNames,
       call=call)
}
#' @rdname tn
#' @aliases tn.tn
#' @method tn tn
#' @export
#' 
#' @examples
## binary vector
tn(c(1, 0, 1, 0, 1))
## Surv object
df0 <- data.frame(t=c(1, 1, 2, 3, 5, 8, 13, 21),
                  e=rep(c(0, 1), 4))
s1 <- with(df0, Surv(t, e, type="right"))
tn(s1, shape="long")
## some awkward values
suppressWarnings(
    s1 <- Surv(time=c(Inf, -1, NaN, NA, 10, 12),
               event=c(c(NA, 1, 1, NaN, Inf, 0.75))))
tn(s1, shape="long", tweOnly=FALSE)
## coxph object
## K&M 2nd ed. Section 1.2. Table 1.1, page 2.
data("hodg", package="KMsurv")
hodg <- data.table::data.table(hodg)
data.table::setnames(hodg,
                     c(names(hodg)[!names(hodg) %in%
                                   c("score", "wtime")],
                       "Z1", "Z2"))
c1 <- coxph(Surv(time=time, event=delta) ~ Z1 + Z2,
            data=hodg[gtype==1 && dtype==1, ])
tn(c1, shape="long")
tn(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
## K&M 2nd ed. Example 7.2, pg 210.
data(kidney, package="KMsurv")
s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
tn(s1)
tn(s1, shape="long")
tn(s1, shape="wide", tweOnly=FALSE)
## formula object
## K&M 2nd ed. Example 7.9, pg 224.
data(kidney, package="KMsurv")
with(kidney, tn(Surv(time=time, event=delta) ~ type, shape="long", twe=F))
## null model
## this is passed to tn.Surv
(t1 <- with(kidney, tn(Surv(time=time, event=delta) ~ 0)))
## but the original call is preserved
attr(t1, "call")
## survival::survfit doesn't accept interaction terms
\dontrun{
    s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)}
## but tn.formula does...
tn(Surv(time=t2, event=d3) ~ z3*z10, data=bmt, shape="wide")
## the same is true for the '.' (dot operator) in formulas
suppressMessages(tn(Surv(time=t2, event=d3) ~ ., data=bmt))
## example where each list element has only one row
## also names are impractical
(t1 <- tn(Surv(time=t2, event=d3) ~ ., data=bmt, sh="long"))
## T&G. Section 3.2, pg 47.
## stratified model
c1 <- coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc)
tn(c1, sh="long")
tn(c1, shape="long")
data(bmt, package="KMsurv")

tn.tn <- function(x, ...,
                  shape=attr(x, "shape"),
                  tweOnly=NULL,
                  abbNames=NULL,
                  call=NULL){
    partMatch(env1=environment(), ...)
    if(shape=="long"){
        x[, "e1" := sum(status), by=list(t, cg)]
        x[, "n1" := max(ncg), by=list(t, cg)]
### x[, "ncg1" := max(ncg), by=list(time, cg)]
        #x[, "status" := NULL]
###
        t1 <- x[, sort(unique(t))]
        lt1 <- length(t1)
        cgInt1 <- x[, unique(as.integer(cg))]
        lcg1 <- length(cgInt1)
        m1 <- '\nMay be inefficient use of memory. Suggest shape="long" instead.\n'
        m2 <- "\nOne covariate group for each time point!"
        if(lcg1 == lt1) message(c(m2, m1))
        m3 <- "\nNote large no. covariate groups relative to no. of time points"
        if(lt1 / lcg1 < 2) message(c(m3, m1))
        ## new memory assignment here
        res1 <- data.table::data.table(
            matrix(rep(as.integer(c(NA, 0)), lcg1 * lt1),
                   ncol=2 * lcg1,
                   nrow=lt1,
                   byrow=TRUE))
        ## import zoo::na.locf.default
        locf <- zoo::na.locf.default
        ## if cg not abbreviated
        ## use as.integer on factor(cg) instead
        abbFn <- if(abbNames){
            identity
        } else {
            as.integer
        }
        ##
        ## a 'for' loop is easier
        ## to read/ debug here
        for(k in seq_along(cgInt1)){
            k1 <- cgInt1[k]
            ## subset by covariate group
            tn1 <- x[abbFn(cg)==k1,
                     list(n1[1L], e1[1]),
                     by="t"]
            ## index of time
            ind1 <- which(t1 %in% tn1[, t])
            ## number at risk
            ## map e.g. k=2 to j=2,3 with
            ## j=2L * k
            j1 <- 2L * k
            data.table::set(res1,
                            i=ind1, j=j1 - 1L,
                            value=tn1[, V1])
            if (ind1[1] > 1L){
                ## set first rows to max(n)
                data.table::set(res1,
                                i=seq.int(ind1[1]), j=j1 - 1L,
                                value=tn1[, max(V1)])
            }
            ## index of missing n
            miss1 <- any(is.na(res1[[j1 - 1L]]))
            if(miss1){
                data.table::set(res1, j=j1 - 1L,
                                value=as.integer(locf(res1[[j1 - 1L]])))
            }
            ## add no. events
            data.table::set(res1,
                            i=ind1, j=j1,
                            value=as.integer(tn1[, V2]))
        }
        ## covariate group names
        if(abbNames){
            cgn1 <- cgInt1
        } else {
            cgn1 <- attr(x, "longNames")[, longName]
        }
        ne_ <- c("n_", "e_")
        ## names for 'n' and 'e' columns
        nne1 <- as.vector(outer(ne_, cgn1, paste, sep=""))
        data.table::setnames(res1, nne1)
        ## make no. at risk (total) per time period
        ## add columns 1, 3, ... , ncol(res1)-1
        res1[, "n" := rowSums(.SD),
             .SDcols = seq.int(from=1L, to=(2L * lcg1 - 1L), by=2L)]
        ## total events per time period
        ## add columns 2, 4, ... , ncol(res1)
        res1[, "e" := rowSums(.SD),
             .SDcols = seq.int(from=2L, to=(2L * lcg1), by=2L)]
        ## now add time
        res1[, "t" := t1]
        data.table::setcolorder(res1,
                                c("t", "n", "e", nne1))
        ln1 <- attr(x, "longNames")
        if(tweOnly){
            ## times with greater than zero events
            x[, "tGr0e" := sum(status) > 0L, by=t]
            ## subset longNames
            cg1 <- x[x[, tGr0e], unique(as.integer(cg))]
            ln1[, "twe" := id %in% cg1]
            ## subset x and drop tGr0e from x
            x <- x[(tGr0e)]
            x[, tGr0e := NULL]
        }
        setAttrTn(res1,
                  "shape"="wide",
                  "tweOnly"=tweOnly,
                  "abbNames"=abbNames,
                  "longNames"=ln1,
                  "ncg"=lcg1,
                  "call"=call,
                  "class"=c("tn", class(x)))
        return(res1)
    }
    if(shape=="wide"){
        n_ <- grep("n_", names(x))
        e_ <- grep("e_", names(x))
        nMe1 <- x[, .SD, , .SDcols=n_] - x[, .SD, ,.SDcols=e_]
        ## no. censored at each time
        e_ <- grep("e_", names(x), value=TRUE)
        substr(e_, 1, 1) <- "c"
        x[, (e_) :=
          nMe1 - data.table::rbindlist(
              list(
                  x[seq.int(2, nrow(x)), .SD, ,.SDcols=n_],
                  as.list(rep(0L, length(n_)))))]
        ## total no. censored per time period
        x[, "nc" := rowSums(.SD),
          .SDcols = grep("c_", names(x))]
        ## no. observations
        nObs1 <- x[,  sum(e + nc)]
        res1 <- data.table::data.table(
            "t" = rep(x[, t], x[, e + nc]))
        res1[, "t" := rep(x[, t], x[, e + nc])]
        res1[, "n" := rev(seq.int(nrow(res1)))]
        res1[, "status" := unlist(mapply(
                            rep, x=c(1,0), times=x[, c(e, nc), by=t]$V1))]
        ## names for covariate groups
        cg_ <- gsub("c_", "", e_)
        res1[, "cg" := unlist(mapply(
                        rep, x=cg_, times=x[, c(e, nc), by=t]$V1))]
        res1[, "cg" := as.factor(cg)]
        ## number at risk per covariate group
        res1[, "ncg" := rev(seq.int(length(n))), by=cg]
        setAttrTn(res1,
                  "shape"="long",
                  "tweOnly"=tweOnly,
                  "abbNames"=abbNames,
                  "longNames"=attr(x, "longNames"),
                  "ncg"=length(cg_),
                  "call"=call,
                  "class"=c("tn", class(x)))
        return(res1)
    }
}
### helper functions
###
## partial matching with an ellipsis
## from environment env1
partMatch <- function(env1=NULL, ...){
    l1 <- as.list(substitute(list(...)))[-1L]
    n1 <- c("sh", "twe", "abb", "con")
    s1 <- sapply(n1, pmatch, names(l1))
    n2 <- c("shape", "tweOnly", "abbNames", "contrasts.arg")
    names(s1) <- n2
    s1 <- s1[!is.na(s1)]
    for (i in seq_along(s1)){
        names(l1)[s1[i]] <- names(s1[i])
    }
    l1 <- l1[names(l1) %in% n2]
    for(i in seq_along(l1)){
        if (is.character(l1[[i]])){
            p1 <- paste0("env1$", names(l1)[i], " <- \"", l1[[i]], "\"")
        } else { 
            p1 <- paste0("env1$", names(l1)[i], " <- ", l1[[i]])
        }
        ## this isn't v. pretty...
        eval(parse(text=p1))
    }
}
## collapse/ paste a data table
## x = data.table
## except = columns to remain unmodified
## nName = new name for collapsed column
## returns the modified data.table
collapseDT <- function(x,
                       except=c("time", "status"),
                       nName="cg"){
### stopifnot(inherits, x, "data.table")
    if(ncol(x)==1){
        data.table::setnames(x, nName)
        return(invisible())
    }
    ## names in 'except'?
    toCollapse1 <- names(x)[!names(x) %in% except]
    x[, (nName) := paste(toCollapse1,
                    .SD,
                    sep="=",
                    collapse=", "),
      .SDcols=toCollapse1,
      by=seq.int(nrow(x))]
    toRemove1 <- which(names(x) %in% toCollapse1)
    if(length(toRemove1)){
        data.table::set(x, j=toRemove1, value=NULL)
    }
    return(invisible())
}
## set attributes for a tn data.table
setAttrTn <- function(x,
                      shape=NULL,
                      tweOnly=NULL,
                      abbNames=NULL,
                      longNames=NULL,
                      ncg=NULL,
                      call=NULL,
                      class=NULL){
### stopifnot(inherits(x, "tn"))
    ## can't use .Internal in a package...
    ## l1 <- .Internal(ls(envir=environment(), all.names=TRUE))
    l1 <- ls()
    l1 <- l1[!grepl("x", l1)]
    for(i in seq_along(l1)){
        data.table::setattr(x,
                            name=l1[i],
                            value=eval(as.name(l1[i])))
    }
    return(x)
}
        toRemove1 <- seq.int(ncol(x))[!seq.int(ncol(x))]
        data.table::set(x, j=xNS1, value=NULL)
        ## collapse strata
        collapseDT(x,
                   except=NA,
                   nName="strat")
   #    
