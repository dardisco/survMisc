#' @name predict
#' @title predicted events
#' @rdname predict
#' @export
#'
predict <- stats::predict
#' @rdname predict
#' @aliases predict.tne
#' @method predict tne
#' @export
#'
#' @param object An objects which has a \code{class} found among:
#'  \cr
#'  methods("predict")
#'  \cr
#' for which an S3 method is available.
#' @param ... Additional arguments (not implemented).
#' @param eMP Add column(s) indicating
#'  \bold{e}vents \bold{m}inus \bold{p}redicted.
#'
#' @details
#' With \eqn{K} covariate groups, We use \eqn{ncg_{ik}}{ncg[i, k]},
#' the number at risk for group \eqn{k},
#' to calculate the number of expected events:
#'  \deqn{P_{ik} = \frac{e_i(ncg_{ik})}{n_i} \quad k=1, 2 \ldots K}{
#'        P[i, k] = e[i] * ncg[i, k] / n[i]}
#'
#' @return The following columns are added to \code{object},
#' by covariate group:
#' \item{P}{predicted number of events}
#' and if \code{eMP==TRUE} (the default)
#' \item{e_P}{observed events minus predicted events}
#' If the \code{tne} object has one column per covariate group
#' (i.e. \code{attr(object, "byWhat")=="time"}), new columns
#' will be
#'
#' @seealso
#' ?survival::predict.coxph
#'
#' @examples
#' ## K&M 2nd ed. Example 7.2, Table 7.2, pp 209--210.
#' data("kidney", package="KMsurv")
#' t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, events.Only=FALSE, by.What="time")
#' predict(t1)
#' stopifnot(t1[, sum(eMP_1 + eMP_2)] <= .Machine$double.neg.eps)
#'
t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, events.Only=FALSE, by.What="time")
c1 <- coxph(Surv(time=time, event=delta) ~ Z1 + Z2,
            data=hodg[gtype==1 && dtype==1, ])
t1 <- tne(c1, by.What="time")

predict(t1)
sf(t1)
#'
predict.tne <- function(object, ..., eMP=TRUE){
    if(attr(object, "by.What")=="status"){
        object[, "P" := (status * ncg) / n]
        ## make events - expected
        if(eMP) object[, "eMP" := status - P]
        return(object)
    }
    if(attr(object, "by.What")=="time"){
        ## names of columns
        n1 <- P1 <- e1 <- grep("n_", names(object), value=TRUE)
        substr(P1, 1L, 2L) <- "P_"
        substr(e1, 1L, 2L) <- "e_"
        eMP1 <- gsub("n_", "eMP_", n1)
        for(i in seq.int(attr(object, "ncg"))){
            object[, (P1)[i] := eval(as.name(n1[i])) * e / n]
            if(eMP){
                object[, (eMP1)[i] :=
                       eval(as.name(e1[i])) - eval(as.name(P1[i]))]
            }
        }
        ## get names for new column order
        tne1 <- c("t", "n", "e")
        ne1 <- c(n1, e1)
        new1 <- c(P1, eMP1[eMP])
        old1 <- names(object)[!names(object) %in%
                              c(tne1, ne1, new1)]
        ## co1 = column order
        co1 <- c(ne1, new1, old1)
        ## new column order
        co1 <- c(tne1, as.vector(t(matrix(co1, nrow=length(n1)))))
        data.table::setcolorder(object, co1)
        return(object)
    }
}
library(data.table)
?setcolorder
?colnames.data.table
