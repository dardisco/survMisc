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
#' t1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, eventsOnly=FALSE, shape="long")
#' predict(t1)
#' stopifnot(attr(t1, "pred")[, sum(eMP_1 + eMP_2)] <= .Machine$double.neg.eps)
#'
#'
predict.tne <- function(object, ..., eMP=TRUE){
    if (attr(object, "shape")=="long"){
        res1 <- matrix(object[, (e * ncg) / n])
        colnames(res1) <- "P"
        ## make events - expected
        if (eMP){
            res1 <- cbind(res1, x[, e - res1])
            colnames(res1) <- c("P", "eMP")
        }
    }
    if(attr(object, "shape")=="wide"){
        ## names of columns
        n1 <- object[, as.matrix(.SD), .SDcols=grep("n_", names(object))]
        res1 <- data.table::data.table(n1 * object[, e / n])
        data.table::setnames(res1, paste("P", seq.int(ncol(res1)), sep="_"))
        if (eMP){
            e1 <- object[, as.matrix(.SD), .SDcols=grep("e_", names(object))]
            na1 <- paste("eMP", seq.int(ncol(res1)), sep="_")
            res1[, (na1) := data.frame(as.matrix(res1) - e1)]
        }}
    data.table::setattr(object, "pred", res1)
    return(attr(object, "pred"))
}
