##' @name tableRhs
##' @export tableRhs
##' @title Table the outcome against all predictors in a formula
##' 
##' @include genSurv.R
##' 
##' @param formula A formula.
##' \cr
##' Works with formulas where the left-hand side is a \code{Surv}
##' object describing right-censored data.
##' @param data A \code{data.frame}.
##' @param return See \bold{Value} below.
##' @param nlf Number of levels defining a factor.
##' \cr
##' Predictors with
##' \eqn{>nlf} levels are considered continuous and are not tabulated.
##' \cr
##' Needs to be less than the number of observations (rows) in the model
##' specified by the \code{formula}.
##' 
##' @return \itemize{
##'  \item If \code{return="summary"} (the default), a \code{table} with
##'  one row per predictor and three columns:
##' \describe{
##'   \item{zeros}{at least one zero present}
##'   \item{someEq}{outcomes equal for least \emph{some} levels of the predictor}
##'   \item{allEq}{outcomes equal for \emph{all} levels of the predictor}
##'   }
##' }
##' Other values return a \code{list} of \code{table}s. Each element is
##' named after the predictor.
##' \itemize{
##' 
##'  \item If \code{return="zeros"}, one \code{table} for each predictor
##'        with a least one zero present.
##'        Each \code{table} shows only those levels
##'        of the predictor for which one level of the outcome is zero.
##' 
##'  \item If \code{return="zEq"}, one \code{table} for each predictor
##'        with a least one \bold{z}ero present or one level which has \bold{eq}ual outcomes.
##'        Each \code{table} shows only those levels where one of the above apply.
##' 
##'  \item If \code{return="counts"}, each \code{table} gives the total
##'        number of levels where zeros and equal outcomes are present and absent.
##' 
##'  \item If \code{return="all"}, a list of \code{table}s of
##'        outcomes for \emph{all} levels of each predictor.
##'
##' }
##' 
##' @details Cross-tabulation of outcomes against levels of a predictor.
##' \cr
##' This is a useful step prior to fitting survival models
##' where the outcome has limited values.
##' 
##' @examples
##' \dontrun{
##' set.seed(1)
##' d1 <- genSurvDf(c=3, rc=0.5, model=FALSE)
##' tableRhs(Surv(t1, e) ~ ., data=d1, return="summary", nlf=2)
##' t1 <- tableRhs(Surv(t1, e) ~ ., data=d1, return="c", nlf=99)
##' ### simple graph
##' p <- par()
##' par( mfrow=c(2, 2))
##' for (i in 1:length(t1)){
##'     graphics::mosaicplot(t1[[i]], main="", cex=1.5)
##'     title(main=list(names(t1[i]), cex=3))
##' }
##' par <- p
##' set.seed(2)
##' d1 <- genSurvDf(f=1, n=10, model=FALSE)
##' t1 <- tableRhs(Surv(t1, e) ~ x1, nlf=9, data=d1)
##' tableRhs(e ~ x1, nlf=9, r="zEq", data=d1)
##' tableRhs(e ~ ., nlf=3, r="c", data=d1)
##' }
tableRhs <- function(formula = y ~ . ,
                     data=parent.frame(),
                     return=c("summary", "zeros", "zEq", "counts", "all"),
                     nlf=2){
### this part taken from stats:lm
    return <- match.arg(return)
    mf <- match.call()
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
### 
    if(nrow(mf) == 0) return("No observations to return")
    if(nlf > nrow(mf)) return("No. levels specified for factor is > no. observations")
### model terms
    mt <- attr(mf, "terms")
### no intercept required for tables
    attr(mt, "intercept") <- 0
### response
    y <- stats::model.response(mf, "numeric")
### is LHS is Surv object for R-censored data
### then split into components
### also get names of outcome variable in Surv object
    if (survival::is.Surv(y)) {
        stopifnot(attr(y, "type")=="right")
        yT <- unname(y[, "time", drop=TRUE])
        yTname <- all.vars(mt)[1]
        y <- unname(y[, "status", drop=TRUE])
        yName<- all.vars(mt)[2]
    } else {
        yName <- all.vars(mt)[1]
    }
### get model matrix
    mm1 <- model.matrix(mt, mf, contrasts)
### keep only columns which are factors
    mm1 <- mm1[, apply(mm1, 2, function(x)
                     length(unique(x)) <= nlf), drop=FALSE]
### check if any columns left
    if (length(mm1)==0) stop("No elements in model.matrix which are factors")
### need SIMPLIFY=FALSE to prevent mapply returning matrix
### i.e. preserve table as class
    l1 <- mapply(table, data.frame(mm1),
                 MoreArgs=list(y, dnn=c("x", yName)),
                 SIMPLIFY=FALSE)
    if (return=="all") return(l1)
### check any with zeros
    z1 <- unlist(lapply(l1, function(x) 0 %in% x ))
### check any duplicates in rows
    eqS <- unlist(lapply(l1, function(x) any(apply(x, 1, anyDuplicated))))
    eqA <- unlist(lapply(l1, function(x) all(x==x[1])))
### 
###----------------------------------------
### summary
###----------------------------------------
### 
    if(return=="summary"){
        n1 <- as.table(matrix(c(z1, eqS, eqA),
                              ncol=length(l1), byrow=TRUE,
                              dimnames=list(c("zeros", "someEq", "allEq"),
                                  names(l1)) ))
        return(n1)
    }
### drop list elments without zeros or some equal
    l1 <- l1[z1 | eqS]
    l1 <- lapply(l1, function(x) {
### ### all elements equal
        eqS <- as.logical(apply(x, 1, anyDuplicated))
        eqA <- all(x==x[1])
### or can use this: eq1 <- apply(x, function(y) all(y==y[1]))
### ### any elements zero
        z1 <- apply(x, 1, function(y) any(y==0))
### ### add columns indicating status
        x <- as.table(cbind(x, "zeros"=z1, "someEq"=eqS))
### ### add table names
        names (dimnames(x)) <- c("x", yName)
        if (eqA) names (dimnames(x)) <- c("x", paste(yName, "all Equal"))
### ### keep only rows with at least one equal or zero
        x <- x[(as.logical(x[, "someEq"] + x[, "zeros"])), ,drop=FALSE]
        return(x)
    })
    if(length(l1)==0) l1 <- NULL
    if(return=="zEq") return(l1)
    if(return=="zeros") {
        l1 <- lapply(l1, function(x) x[x[, "zeros"]==1,
                                       -c((ncol(x)-1), ncol(x))] )
        return(l1)
    }
### get column totals for equal and zeros
    l1 <- lapply(l1, function(x) colSums (x[, (ncol(x)-1):ncol(x), drop=FALSE]))
### double no.s of equal
    #l1 <- lapply(l1, function(x) c(x[1], x[2]*2))
    l1 <- lapply(l1, function(x) cbind("pres"=x, "abs"=nrow(mm1) - x))
###
###----------------------------------------
### else return="counts"
###----------------------------------------
###
    if(length(l1)==0) l1 <- NULL
    return(l1)
}
