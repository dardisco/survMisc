#' @name xtable
#' @title xtable methods
#' @rdname xtable
#' @export
#'
xtable <- function(x,
                   caption=NULL,
                   label=NULL,
                   align=NULL,
                   digits=NULL,
                   display=NULL,
                   ...){
    UseMethod("xtable")
}
#' @rdname xtable
#' @aliases xtable.Surv
#' @method xtable Surv
#' @export
#'
#' @param x An objects which has a \code{class} found among
#'  \cr
#'  methods("xtable")
#'  \cr
#' for which an S3 method is available.
#' @param ... Additional arguments
#'
#' @details
#' @return An object of \code{class} "xtable",
#'  which inherits the \code{data.frame}
#'  class and contains several additional attributes specifying the
#'  table formatting options.#'
#'
#' @examples
#' ## K&M 2nd ed. Example 7.2, Table 7.2, pp 209--210.
#' data("kidney", package="KMsurv")
#' t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney)
#' E(t1)
#' stopifnot(t1[, sum(eME_1 + eME_2)] <= .Machine$double.neg.eps)
#'
set.seed(1)
s1 <- Surv(time=rep(1:8, 8),
           event=sample(c(0,1), size=64, replace=TRUE))
xtable(s1)
xtable.Surv <- function(x,
                        caption=paste0("\\texttt{Surv} object, ",
                        attr(x, "type"),
                        " censored"),
                        label=NULL,
                        align=NULL,
                        digits=NULL,
                        display=NULL,
                        ...){
    ## from as.character.default()
    x <- .Internal(as.vector(x, "character"))
    xtable(x,
           caption=caption,
           label=label,
           align=align,
           digits=digits,
           display=display,
           ...)
}
options(xtable.width = NULL)
options(xtable.counter = 1)
c1 <- as.character(s1)
xtable.character <- function(x,
                             caption=NULL,
                             label=TRUE,
                             align=c("l", rep("c", nCol)),
                             digits=NULL,
                             display=NULL,
                             ...,
                             nCol=getOption("survMisc.nCol", 10)){
    m1 <- matrix(c(x, rep("", length(x) %/% nCol)),
                 ncol=nCol,
                 byrow=TRUE)
    if(is.logical(label) && label) label <- labelF()
    x1 <- xtable(m1,
                 caption=caption,
                 label=label,
                 align=align,
                 digits=digits,
                 display=display)
    print(x1,
          include.colnames=getOption("xtable.include.colnames",
          FALSE),
          booktabs=getOption("xtable.booktabs", TRUE),
          ...)
    return(invisible(x1))
}
set.seed(1)
t1 <- tne(sample(c(0,1), 200, replace=TRUE))
xtable(t1)
data("pbc", package="survival")
t1 <- tne(Surv(time, status==2) ~ age, data=pbc, byWhat="time")
xtable.data.table <- function(x,
                              caption=paste0("\\texttt{", attr(x, "call"), "}"),
                              label=TRUE,
                              align=TRUE,
                              digits=NULL,
                              display=NULL,
                              ...){
    if(nrow(x) > getOption("datatable.print.nrows")){
        pt1 <- getOption("datatable.print.topn")
        x[, "row" := seq.int(nrow(x))]
        data.table::setcolorder(x, c("row", names(x)[1:ncol(x)-1L]))
        spc1 <- as.list(c("---", rep("", (ncol(x) - 1))))
        d1 <- data.table::rbindlist(
            list(x[seq.int(pt1), ],
                 spc1,
                 x[rev((nrow(x) + 1L) - seq.int(pt1)), ]))
    } else {
        d1 <- x
    }
    if(is.logical(label) && label) label <- labelF()
    data.table::setattr(d1, "class", "data.frame")
    x1 <- xtable(d1,
                 caption=caption,
                 label=label,
                 align=NULL,
                 digits=digits,
                 display=display
                 ...)
    print.xtable(x1,
                 include.rownames=getOption("xtable.include.rownames",
                 FALSE),
                 ...)
    return(invisible(x1))
}
s1 <- survfit(Surv(time, status==2) ~ strata(edema) + edema, data=pbc)
survfit(c1)
xtable.survfit <- function(x,
                           caption=paste0("\\texttt{", x$call, "}")
                           label=TRUE,
                           align=TRUE,
                           digits=NULL,
                           display=NULL,
                           ...){
### based on
### survival:::survmean
    sm <- survival:::survmean
    debugonce("sm")
    sm(x, rmean=1000)
    colN1 <- c("records", "n.start", "events",
               "median", paste(x$conf.int, c("LCL", "UCL"),
                               sep = ""))
    rowN1 <- names(x$strata)
}

###----------------------------------------
### helper functions
labelF <- function(){
    if(is.null(getOption("xtable.counter"))){
        options(xtable.counter = 1)
    } else {
        options(xtable.counter = getOption("xtable.counter") + 1)
    }
    return(paste0("tab:", getOption("xtable.counter")))
}
