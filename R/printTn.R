#' @title 'print' methods
#' @rdname printTn
#' @details Prints a \code{tn} object with 'nice' formatting.
#'  \cr
#'  Options may be set for a session using e.g.
#'  \cr
#'  options(survMisc.nColSP=4L)
#'  \cr
#'  It is similar to the behavior of \code{print.data.table} but
#'  has additional arguments controlling the number of columns
#'  sent to the terminal.
#' @param x An object of class \code{tn}.
#' @param ... Additional arguments (not implemented).
#' @param maxRow Maximum number of rows to print.
#'  \cr
#'  If \code{nrow(x) > maxRow}, just the first and last
#'  \code{nRowP} are printed.
#' @param nRowP \bold{N}umber of rows to \bold{p}rint from
#'  the start and end of the object.
#' @param pRowNames Print row names? Default is \code{TRUE}.
#' @param maxCol Maximum number of columns to print.
#'  \cr
#'  If \code{ncol(x) > maxCol}, just the first \code{nColSP}
#'  and last \code{maxCol - nColSP} columns are printed.
#' @param nColSP \bold{N}umber of \bold{col}umns to \bold{p}rint from
#'  the \bold{s}tart of the object.
#' @param sigDig \bold{Sig}nificant \bold{dig}its.
#'  \cr
#'  This is passed as an argument to
#'  \cr
#'  ?signif
#'  when preparing the object for printing.
#' @return A printed representation of the \code{tn} object
#'  is send to the terminal as a \emph{side effect} of
#'  calling the function.
#'  \cr
#'  The return value cannot be \code{assign}ed.
#' @author Based on existing work by Brian Diggs.
#' @seealso
#' data.table:::print.data.table
#' @note
#' All numeric arguments to the function must be supplied as integers.
#'
#' @aliases print.tn
#' @method print tn
#' @export
#'
#' @examples
#' set.seed(1)
#' (x <- data.table::data.table(matrix(rnorm(1800), ncol=15, nrow=120)))
#' data.table::setattr(x, "class", c("tn", class(x)))
#' p1 <- print(x)
#' stopifnot(is.null(p1))
#' x[1:80, ]
#' x[0, ]
#' (data.table::set(x, j=seq.int(ncol(x)), value=NULL))
print.tn <- function(x, ...,
                      maxRow=getOption("datatable.print.nrows", 100L),
                      nRowP=getOption("datatable.print.topn", 5L),
                      pRowNames=TRUE,
                      maxCol=getOption("survMisc.maxCol", 8L),
                      nColSP=getOption("survMisc.nColSP", 7L),
                      sigDig=getOption("survMisc.sigDig", 2L)){
    if(nrow(x)==0L){
        if (length(x)==0L){
            cat("Null 'tn' object (0 rows and 0 cols)\n")
        } else {
            cat("Empty 'tn' object (0 rows) of ", length(x),
                " col", if (length(x) > 1L)
                "s", ": ", paste(head(names(x), 6), collapse = ", "),
                if (ncol(x) > 6)
                "...", "\n", sep = "")
        }
        return(invisible())
    }
    stopifnot(all(sapply(list(maxRow, nRowP, maxCol, nColSP),
                         as.integer)))
    stopifnot(maxRow > nRowP)
    ## lCol1 = last columns; needs to be at least one
    stopifnot((lCol1 <- maxCol - nColSP) >= 1)
    if(nrow(x) > maxRow){
        toPrint1 <- rbind(head(x, nRowP),
                          tail(x, nRowP))
        ## row names
        rn1 <- c(seq_len(nRowP),
                 seq.int(to=nrow(x), length.out=nRowP))
        rowDots1 <- TRUE
    } else {
        toPrint1 <- x
        rn1 <- seq_len(nrow(x))
        rowDots1 <- FALSE
    }
    if(ncol(x) > (nColSP + lCol1 + 1L)){
        toPrint1 <- cbind(
            toPrint1[, seq.int(nColSP), with=FALSE],
            toPrint1[, seq.int(to=ncol(x), length.out=lCol1), with=FALSE])
        colDots1 <- TRUE
    } else {
        colDots1 <- FALSE
    }
    toPrint1 <- do.call("cbind",
                        lapply(toPrint1,
                               function(col) signif(col,
                                                    digits=sigDig)))
    if(pRowNames){
        rownames(toPrint1) <- paste(format(rn1, right = TRUE),
                                    ":",
                                    sep = "")
    } else {
        rownames(toPrint1) <- rep.int("", nrow(x))
    }
    if(rowDots1){
        toPrint1 <- rbind(head(toPrint1, nRowP),
                          "---" = "",
                           tail(toPrint1, nRowP))
        rownames(toPrint1) <- format(rownames(toPrint1),
                                     justify="right")
    }
    if(colDots1){
        toPrint1 <- cbind(
            toPrint1[, seq.int(nColSP), drop=FALSE],
            rep("", nrow(toPrint1)),
            toPrint1[, seq.int(to=ncol(toPrint1), length.out=lCol1),
                     drop=FALSE])
        colnames(toPrint1)[colnames(toPrint1)==""] <- " ---"
    }
    if(!rowDots1){
        toPrint1 <- rbind(toPrint1,
                          matrix(colnames(toPrint1), nrow=1L))
    }
    print(toPrint1, right=TRUE, quote=FALSE)
    return(invisible())
}
