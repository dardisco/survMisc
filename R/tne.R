#' @name tne
#' @title \bold{t}ime and \bold{n}umber at risk, with status or
#' numer of events.
#'
#' @include print.R
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
#' The method for \code{tne} objects converts the object
#' from 'long' to 'wide' form or \emph{vice versa}.
#' @param shape See \bold{Value} below.
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
#' as an \code{attribute} of the returned \code{tne} object.
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
#'  to the final \code{tne.data.table} method. 
#'
#' @return A \code{data.table} with the additional \code{class}
#'  \code{tne}. The same information can be displayed in two ways.
#'  \cr
#' If \code{shape=="wide"} (the default), the data is
#'  returned in 'wide' format.
#'  \cr
#' There is one row for each time.
#'  \cr
#' For a \code{numeric} or \code{Surv} object this has columns: 
#'  \item{t}{time.}
#'  \item{n}{number at risk.}
#'  \item{e}{Number of events.}
#'  A \code{survfit}, \code{coxph} or \code{formula} object
#'  will have additional columns for \code{n} and \code{e}
#'  for each covariate group. 
#'
#' If \code{shape=="long"}, the data is
#'  returned in 'long' format.
#'  \cr
#' There is one row for each unique timepoint per covariate group.
#'  \cr
#' For a \code{numeric} or \code{Surv} object this has columns: 
#'  \item{time}{time}
#'  \item{n}{number at risk}
#'  \item{e}{}
#'  A \code{survfit}, \code{coxph} or \code{formula} object
#'  will have additional columns:
#'   \item{cg}{\bold{c}ovariate \bold{g}roup.
#'  \cr
#'  This is formed from combining the variables; these
#'  are separated by an comma ','}
#'   \item{ncg}{\bold{n}umber at risk, by \bold{c}ovariate \bold{g}roup}
#'  
#' \bold{Special terms}.
#'  \cr \cr
#' The following are considered 'special'
#' terms in a survival model:
#'   \item{strata}{For a stratified model, a \code{by list} with
#'    one \code{tne} object per strata as above.}
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
#' @rdname tne
#' @export
#'
## for R CMD check
n <- e <- status <- NULL
tne <- function(x, ...) UseMethod("tne")
## all are methods ultimately passed to
## tne.data.frame (below)
## except tne.numeric() and tne.Surv()
###----------------------------------------
#' @rdname tne
#' @aliases tne.numeric
#' @method tne numeric
#' @export
#'
tne.numeric <- function(x, ...){
    stopifnot(all(x >= 0 && x <=1))
    res1 <- data.table::data.table(
        "t"=(t <- seq_along(x)),
        "n"=rev(t),
        "e"=x)
    setAttrTne(res1,
               shape="long",
               abbNames=TRUE,
               ncg=0,
               call=match.call(),
               class=c("tne", class(res1)))
    return(res1)
}
#' @rdname tne
#' @aliases tne.Surv
#' @method tne Surv
#' @export
#'
tne.Surv <- function(x, ..., 
                    call=NULL){
    stopifnot(inherits(x, "Surv"))
    stopifnot(attributes(x)$type=="right")
    if(is.null(call)) call <- match.call()
    res1 <- data.table::data.table(unclass(x))
    data.table::setkey(res1, "time")
    res1 <- res1[, list("n"=length(status),
                        "e"=sum(status)),
                 by=sort(time, na.last=TRUE)]
    res1[, "n" := c(sum(n), sum(n) - cumsum(n)[ - length(n)])]
    data.table::setnames(res1, c("t", "n", "e"))
    setAttrTne(res1,
               shape="long",
               ncg=0,
               call=call,
               class=c("tne", class(res1)))
    return(res1)
}
###
###----------------------------------------
###
#' @rdname tne
#' @aliases tne.coxph
#' @method tne coxph
#' @export
#'
tne.coxph <- function(x, ...,
                     shape=c("wide", "long"),
                      abbNames=TRUE,
                      contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    x$call$formula <- stats::terms(
        x=formula(x$call),
        specials=c("strata", "cluster", "tt"))
    call1 <- x$call
    x$call$drop.unused.levels <- TRUE
    x$call[[1]] <- as.name("model.frame")
    ## model.frame
    xMF1 <- eval(x$call, parent.frame())
    tne(x=xMF1,
        shape=shape,
        abbNames=abbNames,
        contrasts.arg=contrasts.arg,
        call=call1)
}
#' @rdname tne
#' @aliases tne.survfit
#' @method tne survfit
#' @export
#'
tne.survfit <- function(x, ..., 
                        shape=c("wide", "long"),
                        abbNames=TRUE,
                        contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    x$call$formula <- stats::terms(
        x=formula(x$call),
        specials=c("strata", "cluster", "tt"))
    call1 <- x$call
    x$call$drop.unused.levels <- TRUE
    x$call[[1]] <- as.name("model.frame")
    xMF1 <- eval(x$call, parent.frame())
    tne(x=xMF1,
        shape=shape,
        abbNames=abbNames,
        contrasts.arg=contrasts.arg,
        call=call1)
}
#' @rdname tne
#' @aliases tne.formula
#' @method tne formula
#' @export
#'
tne.formula <- function(x, ...,
                        shape=c("wide", "long"),
                        abbNames=TRUE,
                        contrasts.arg=NULL){
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    stopifnot(inherits(x, "formula"))
### based on code from stats::lm()
    mc1 <- match.call()
    names(mc1)[names(mc1)=="x"] <- "formula"
    mc1 <- mc1[c(1L, match(c("formula", "data"), names(mc1), 0L))]
    call1 <- mc1
    mc1$drop.unused.levels <- TRUE
    mc1[[1L]] <- as.name("model.frame")
    mf1 <- eval(mc1, parent.frame())
    tne(x=mf1,
        shape=shape,
        abbNames=abbNames,
        contrasts.arg=contrasts.arg,
        call=call1)
}
#' @rdname tne
#' @aliases tne.data.frame
#' @method tne data.frame
#' @export
#'
tne.data.frame <- function(x, ...,
                           shape=c("wide", "long"),
                           abbNames=TRUE,
                           contrasts.arg=NULL,
                           call=NULL){
    stopifnot(survival::is.Surv(x[[1]]))
    stopifnot(attr(x[[1]], "type") == "right")
    shape <- match.arg(shape)
    partMatch(env1=environment(), ...)
    ## data.table from x
    xDT <- data.table::data.table(
        cbind(
            stats::model.matrix(terms(x), x,
                                contrasts.arg=contrasts.arg),
            stats::model.response(x)))
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
        xDTSn1 <- grep("^strata\\(.*\\)", names(xDT))
        ## separate table onlly for strata
        xDTstr <- xDT[, .SD, .SDcols=xDTSn1]
        collapseDT(xDTstr, except=NA, nName="strat")
        data.table::set(xDT, j=xDTSn1, value=NULL)
        collapseDT(xDT,
                   except=c("time", "status"),
                   nName="cg")
        xDT[, "cg" := as.factor(cg)]
        xDT[, "strat" := as.factor(xDTstr[, strat])]
                                        #browser()
        xDTSn1 <- seq.int(names(xDT))[!(grepl("^strat", names(xDT)))]
        res1 <- vector(mode="list")
        for (i in xDT[, levels(strat)]){
            res1[[i]] <- xDT[levels(strat)==i, .SD, .SDcols=xDTSn1]
        }
        lapply(res1, function(x) x[, "c":=1])
        lapply(res1, tne,
               shape=shape,
               abbNames=abbNames)
        names(res1) <- xDT[, levels(strat)]
        data.table::setattr(res1,
                            "class",
                            c("stratTne", class(res1)))
        data.table::setattr(res1, "call", call)
        return(res1)
    } else {
        if(stats::is.empty.model(x)){
            ## convert to Surv object
            s1 <- xDT[, survival::Surv(time, status)]
            return(tne(s1,
                       call=call))
        } else {
            collapseDT(xDT,
                       except=c("time", "status"),
                       nName="cg")
            tne(x=xDT,
                shape=shape,
                abbNames=abbNames,
                call=call)
        }
    }
}
#' @rdname tne
#' @aliases tne.data.table
#' @method tne data.table
#' @export
#' 
tne.data.table <- function(x, ...,
                          shape=c("wide", "long"),
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
        "id" = x[, seq_along(levels(cg))],
        "longName" = x[, levels(cg)])
    if(abbNames) x[, "cg" := as.integer(cg)]
    data.table::setnames(x, c("t", colnames(x)[-1]))
    x[, "e" := sum(status), by=list(t, cg)]
    x[, "ncg" := max(ncg), by=list(t, cg)]
    x[, "n" := max(n), by=list(t)]
    x[, "status" := NULL]
    x <- x[!(duplicated(x)), ]
    data.table::setcolorder(x,
                            c("t", "n", "e", "cg", "ncg"))
    if(shape=="wide"){
        setAttrTne(x,
                   "shape"="long",
                   "abbNames"=abbNames,
                   "longNames"=ln1,
                   "ncg"=nrow(ln1),
                   "call"=call,
                   "class"=c("tne", class(x)))
        return(tne(x=x,
                   shape="long",
                   abbNames=abbNames,
                   call=call))
    } else {
        setAttrTne(x,
                   "shape"="long",
                   "abbNames"=abbNames,
                   "longNames"=ln1,
                   "ncg"=nrow(ln1),
                   "call"=call,
                   "class"=c("tne", class(x)))
        data.table::setkey(x, "cg")
        return(x)
    }
}
#' @rdname tne
#' @aliases tne.tne
#' @method tne tne
#' @export
#' 
#' @examples
##' ## binary vector
##' tne(c(1, 0, 1, 0, 1))
##' ## Surv object
##' df0 <- data.frame(t=c(1, 1, 2, 3, 5, 8, 13, 21),
##'                   e=rep(c(0, 1), 4))
##' s1 <- with(df0, Surv(t, e, type="right"))
##' tne(s1, shape="long")
##' ## some awkward values
##' suppressWarnings(
##'     s1 <- Surv(time=c(Inf, -1, NaN, NA, 10, 12),
##'                event=c(c(NA, 1, 1, NaN, Inf, 0.75))))
##' tne(s1)
##' ## coxph object
##' ## K&M 2nd ed. Section 1.2. Table 1.1, page 2.
##' data("hodg", package="KMsurv")
##' hodg <- data.table::data.table(hodg)
##' data.table::setnames(hodg,
##'                      c(names(hodg)[!names(hodg) %in%
##'                                    c("score", "wtime")],
##'                        "Z1", "Z2"))
##' c1 <- coxph(Surv(time=time, event=delta) ~ Z1 + Z2,
##'             data=hodg[gtype==1 && dtype==1, ])
##' tne(c1, shape="long")
##' tne(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
##' ## K&M 2nd ed. Example 7.2, pg 210.
##' data("kidney", package="KMsurv")
##' with(kidney[kidney$type==2, ], tne(Surv(time=time, event=delta), shape="long"))
##' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' tne(s1)[e > 0, ]
##' ## A null model is passed to tne.Surv
##' (t1 <- with(kidney, tne(Surv(time=time, event=delta) ~ 0)))
##' ## but the original call is preserved
##' attr(t1, "call")
##' ## survival::survfit doesn't accept interaction terms...
##' \dontrun{
##'     s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)}
##' ## but tne.formula does:
##' tne(Surv(time=t2, event=d3) ~ z3*z10, data=bmt, shape="wide")
##' ## the same is true for the '.' (dot operator) in formulas
##' tne(Surv(time=t2, event=d3) ~ ., data=bmt)
##' ## better to store this with shape=="long"
##' ## also names are impractical
##' (t1 <- tne(Surv(time=t2, event=d3) ~ ., data=bmt, sh="long"))
##' attr(t1, "longNames")
##' ## T&G. Section 3.2, pg 47.
##' ## stratified model
##' c1 <- coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc)
##' tne(c1, sh="long")
##' tne(c1, shape="long")
tne.tne <- function(x, ...,
                    shape=attr(x, "shape"),
                    abbNames=NULL,
                    call=NULL){
    partMatch(env1=environment(), ...)
    if (shape=="long"){
        t1 <- x[, sort(unique(t))]
        ## length of time (number of observations)
        lt1 <- length(t1)
        ## covariate groups (as integer)
        cgInt1 <- x[, seq_along(unique(cg))]
        ## length (number) of covariate groups
        lcg1 <- length(cgInt1)
        m1 <- '\nMay be inefficient use of memory. Suggest shape="long" instead.\n'
        m2 <- "\nOne covariate group for each time point!"
        if(lcg1 == lt1) message(c(m2, m1))
        m3 <- "\nNote large no. covariate groups relative to no. of time points"
        if(lt1 / lcg1 < 2) message(c(m3, m1))
### new structure here
        res1 <- data.table::data.table(
            matrix(rep(as.integer(c(NA, 0)), lcg1 * lt1),
                   ncol=2 * lcg1,
                   nrow=lt1,
                   byrow=TRUE))
        ## import zoo::na.locf.default
        locf <- zoo::na.locf.default
        ## abbFn = abbreviate function
        ## if cg not abbreviated
        ## use as.integer on factor(cg) instead
        if (abbNames){
          abbFn <- identity
        } else {
          abbFn <- as.integer
        }
        ## a 'for' loop is easier
        ## to read/ debug here
        for (k in cgInt1){
            ## subset by covariate group
            tne1 <- x[abbFn(cg)==k,
                     list(ncg[1L], e[1]),
                     by="t"]
            ## index of time
            ind1 <- which(t1 %in% tne1[, t])
            ## number at risk
            ## map e.g. k=2 to j=2,3 with
            ## j=2L * k
            j1 <- 2L * k
            data.table::set(res1,
                            i=ind1, j=j1 - 1L,
                            value=tne1[, V1])
            if (ind1[1] > 1L){
                ## set first rows to max(n)
                data.table::set(res1,
                                i=seq.int(ind1[1]), j=j1 - 1L,
                                value=tne1[, max(V1)])
            }
            ## index of missing n
            miss1 <- any(is.na(res1[[j1 - 1L]]))
            if (miss1){
                data.table::set(res1,
                                i=seq.int(max(ind1)),
                                j=j1 - 1L,
                                value=as.integer(locf(res1[[j1 - 1L]], fromLast=TRUE)))
                if (max(ind1) < nrow(res1)){
                    data.table::set(res1,
                                    i=seq.int(from=max(ind1) + 1, to=nrow(res1)),
                                    j=j1 - 1L,
                                    value=0)
                }
            }
            ## add no. events
            data.table::set(res1,
                            i=ind1, j=j1,
                            value=as.integer(tne1[, V2]))
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
        ln1 <- data.table::copy(attr(x, "longNames"))
        setAttrTne(res1,
                  "shape"="wide",
                  "abbNames"=abbNames,
                  "longNames"=attr(x, "longNames"),
                  "ncg"=lcg1,
                  "call"=call,
                  "class"=c("tne", class(res1)))
        return(res1)
    }
    if(shape=="wide"){
        n_ <- grep("n_", names(x))
        e_ <- grep("e_", names(x))
        ## no. at risk - no. events
        nMe1 <- x[, .SD, .SDcols=n_] - x[, .SD, .SDcols=e_]
        ## add no. censored
        c1 <- nMe1[seq.int(nrow(x) - 1), ] - x[seq.int(2, nrow(x)), .SD, .SDcols=n_]
        c1 <- data.table::rbindlist(list(
            c1,
            x[nrow(x), .SD, .SDcols=n_]))
        ## names for censored columns
        c_ <- grep("e_", names(x), value=TRUE)
        substr(c_, 1, 1) <- "c"
        x[, (c_) := c1]
        ## total no. censored per time period
        x[, "nc" := rowSums(.SD),
          .SDcols = grep("c_", names(x))]
        ## names of covariate groups
        n1 <- sub("c_", "", c_)
        l1 <- vector(mode="list", length=length(n1))
        for (i in seq_along(n1)){
            e_n <- grep(paste0("e_", n1[i]), names(x))
            c_n <- grep(paste0("c_", n1[i]), names(x))
            ## which times have at least one event
            ## or one censored observation
            l1[[i]] <- as.logical(rowSums(x[, .SD, .SDcols=c(e_n, c_n)]))
        }
### new structure here
        res1 <- data.table::data.table(
            t=unlist(lapply(l1, function(i) x[i, t])),
            n=unlist(lapply(l1, function(i) x[i, n])),
            e=unlist(sapply(seq_along(n1),
                function(i) x[, .SD, .SDcols=e_[i]][l1[[i]], ])),
            cg=as.integer(unlist(mapply(rep, n1, sapply(l1, sum)))),
            ncg=unlist(sapply(seq_along(n1),
                function(i) x[, .SD, .SDcols=n_[i]][l1[[i]], ])))
        data.table::setkey(res1, t)
        setAttrTne(res1,
                  "shape"="long",
                  "abbNames"=abbNames,
                  "longNames"=attr(x, "longNames"),
                  "ncg"=length(cg_),
                  "call"=call,
                  "class"=c("tne", class(x)))
        return(res1)
    }
}
### helper functions
###
## partial matching with an ellipsis
## from environment env1
partMatch <- function(env1=NULL, ...){
    stopifnot(is.environment(env1))
    l1 <- as.list(substitute(list(...)))[-1L]
    n1 <- c("sh", "abb", "con")
    s1 <- sapply(n1, pmatch, names(l1))
    n2 <- c("shape", "abbNames", "contrasts.arg")
    names(s1) <- n2
    s1 <- s1[!is.na(s1)]
    for (i in seq_along(s1)){
        names(l1)[s1[i]] <- names(s1[i])
    }
    l1 <- l1[names(l1) %in% n2]
    for(i in seq_along(l1)){
        ## this isn't v. pretty...
        if (is.character(l1[[i]])){
            p1 <- paste0("env1$", names(l1)[i], " <- \"", l1[[i]], "\"")
        } else { 
            p1 <- paste0("env1$", names(l1)[i], " <- ", l1[[i]])
        }
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
    stopifnot(inherits(x, "data.table"))
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
## set attributes for a tne data.table
setAttrTne <- function(x,
                       shape=NULL,
                       abbNames=NULL,
                       longNames=NULL,
                       ncg=NULL,
                       call=NULL,
                       class=NULL){
### stopifnot(inherits(x, "tne"))
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
