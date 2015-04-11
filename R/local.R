##' @name local
##' @rdname local
##' @title Local tests for a model
##' @keywords htest
##' 
##' @export locScore
##' 
locScore <- function(x, ...){
    UseMethod("locScore")
}
##' @rdname local
##' @aliases locScore.coxph
##' @method locScore coxph
##' @export
##'
##' @param x A model of class \code{coxph}
##' @param ... Additional arguments (not implemented)
##' @param all Fit \emph{all} combinations of predictors
##' @param hypo Hypothesis to test.
##' \cr
##' There should be at least one coefficient to exclude and one to keep.
##' \cr
##' This is a specified as vector of the same length as the number of
##' coefficients in the model.
##' \cr
##' This should be a logical vector
##' (i.e. composed of \code{TRUE} and \code{FALSE} or a vector of
##' \eqn{0}s and \eqn{1}s.
##' \itemize{
##'   \item \code{FALSE} or zeros indicate coefficients to exclude
##'   \item \code{TRUE} or ones indicate coefficients to keep.
##' }
##' @param ties Method of handling ties when refitting model.
##' \cr
##' Must be one of \code{breslow}, \code{efron} or \code{exact}.
##' @return For \code{locScore} a \code{list} with the following elements,
##' which are \code{data.table}s:
##' \item{coef}{coefficients from refitted model(s)}
##' \item{score}{hypothesis and chi-square test}
##' For \code{locLR} and \code{locWald}, a \code{data.table}
##' showing the hypothesis being tested and the results of the test.
##'
##' @details The null hypothesis is that some of the coefficients in
##' the model are zero (\eqn{H_0 :  \hat{B}_i=0, \quad i \geq 1}{H0: Bi=0})
##' vs. the alternative that at least one of them is nonzero.
##' \cr \cr
##' All of these tests are distributed as chi-square
##' with degrees of freedom \eqn{=} number of excluded coefficients.
##' \cr \cr
##' For the \bold{score} test,
##' the model is fitted again with the coefficients of
##' interest excluded.
##' \cr
##' A value for the remaining coefficients is obtained.
##' Then the complete model is fit again using these new values as
##' initial values for
##' those remaining coefficients and using zero as the initial
##' value for the excluded coefficients.
##' \cr
##' Values for the excluded coefficients are generated without
##' iteration. (I.e. the first values calculated, with no convergence towards
##' maximum likelihood estimators).
##' \cr
##' The test is:
##' \deqn{ \chi_{SC}^2 = U^T I^{-1} U}{
##'  SC = U I^-1 U}
##' where \eqn{U} is the score vector  and \eqn{I^{-1}}{I^-1} is the
##' covariance or inverse of the information matrix.
##' (These are given by \code{colSums(survival::coxph.detail(x)$score)}
##' and \code{x$var} respectively).
##' \cr \cr
##' For the \bold{likelihood ratio} test, the model is also refit with
##' the coefficients of interest excluded. The likelihood ratios
##' from the full model and those with coefficients excluded
##' are used to construct the test:
##' \deqn{ \chi^2_{LR} = 2(LR_{full} - LR_{excluded})}{
##'  LR = 2*( LR[full] - LR[excluded] )}
##' \cr \cr
##' The \bold{Wald} chi-squared statistic is given by:
##' \deqn{ \chi^2_W = \hat{B}^T I^{-1} \hat{B} }{
##'  W = B I^-1 B }
##' Where \eqn{\hat{B}}{B} is the vector of fitted coefficients
##' (from the complete model) thought to be \eqn{=0}.
##' \cr
##' \eqn{I^{-1}}{I^-1} is composed of the corresponding elements from the
##' covariance matrix of the model.
##' @examples
##' data(larynx, package="KMsurv")
##' c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx)
##' locScore(c1, all=TRUE)
##' locScore(c1, hypo=c(0, 0, 0, 1))
##' locScore(coxph(Surv(time, delta) ~ stage + age, data=larynx))
##' @references Examples are from:
##' \bold{K&M} Example 8.2, pp 264-6.
locScore.coxph <- function(x,
                           ...,
                           all=FALSE,
                           hypo=NULL,
                           ties=c("breslow", "efron", "exact")
                           ){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### get names of coefficients
    n1 <- names(x$coefficients)
### predictor terms
    ties <- match.arg(ties)
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
### if no hypothesis, set first coefficient to zero
    if (is.null(hypo)) hypo <- c(0, rep(1, length=(np1-1)))
    if (length(hypo) != np1) stop ("Hypothesis must be same length as no. of coefficients")
    m1 <- "Hypothesis must be a logical vector of TRUE/FALSE (or 0/1) statements"
    if (!all(c(0, 1) %in% hypo | !all(c(TRUE, FALSE) %in% hypo))) stop (m1)
### taken from code for lm()
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
    if (!all){
        r1 <- .getScore(DT=dt1, hypo=hypo,
                        mFrame=mf, x=x, ties=ties)
### r1$score[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
        return(r1)
    }
###
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
### drop first and last rows (i.e. all present or all absent)
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
    colnames(c2) <- n1
    sc1 <- data.table(c2)
    co1 <- data.table(c2)
### add additional columns to store values
### give default value as double
    sc1[, c("chiSq", "df", "pVal") := 0.1]
    for (i in seq(nrow(c2))){
        hypo <- c2[i, ]
        r1 <- .getScore(DT=dt1, hypo=hypo,
                        mFrame=mf, x=x, ties=ties)
        sc1[i, ] <- r1$score
        set(co1, i=i, j=1:ncol(co1), value=r1$coef)
    }
    sc1 <- sc1[order(sc1$pVal), ]
### sc1[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
    return(list(
        coef=co1,
        score=sc1
        ))
}
###
.getScore <- function(DT, hypo, mFrame, x, ties){
### for R CMD check
    chiSq <- NULL
### get names of coefficients
    n1 <- names(x$coefficients)
### copy to modify (drop terms later)
    dt2 <- copy(DT)
    n2 <- n1[as.logical(hypo)]
    cn1 <- seq_along(colnames(dt2))
### get position of columns not in hypothesis
    i1 <- as.integer(cn1[!colnames(dt2) %in% n2])
### drop these unused columns
    set(dt2, j=i1, value=NULL)
    stopifnot(attr(model.response(mFrame), "type")=="right")
    y1 <- data.table(unclass(model.response(mFrame, "numeric")))
    dt2[, c("t", "e") := y1]
    coef1 <- coxph(Surv(t, e) ~ ., ties=ties, data=dt2)$coefficients
### add to results
    init1 <- as.vector(hypo)
    init1[hypo==1] <- coef1
### get original formula
### then refit with new initial values and no iterations
### (i.e. no convergence of estimates)
### make formula look up the eval() tree
### rather than just where it was called originally
    attr(x$formula, ".Environment") <- parent.frame()
    cox1 <- survival::coxph(x$formula,
                            data= eval(x$call$data),
                            ties=ties,
                            init=init1,
                            iter.max=0)
### score vector
    sc1 <- colSums(survival::coxph.detail(cox1)$score)
### results
###  init1 <- formatC(init1)
    r1 <- data.table(t(init1))
    setnames(r1, n1)
    r2 <- data.table(t(hypo),
                     sc1 %*% cox1$var %*% sc1,
                     sum(hypo)
                     )
    setnames(r2, c(n1, "chiSq", "df"))
    r2[, "pVal" := (1-stats::pchisq(chiSq, df))]
###
    return(list(coef=r1,
                score=r2))
}
###----------------------------------------
###
##' @rdname local
##' @aliases locLR
##' @export 
##'
locLR <- function(x, ...){
    UseMethod("locLR")
}
##' @rdname local
##' @aliases locLR.coxph
##' @method locLR coxph
##' @export
##' 
##' @examples
##' ###
##' data(larynx, package="KMsurv")
##' c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx, method="breslow")
##' locLR(c1, all=TRUE)
##' locLR(c1, hypo=c(FALSE, FALSE, FALSE, TRUE))
##'
locLR.coxph <- function(x, ...,
                        all=FALSE,
                        hypo=NULL,
                        ties=c("breslow", "efron", "exact")
                        ){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### get names of coefficients
    n1 <- names(x$coefficients)
    ties <- match.arg(ties)
### no. predictor terms
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
### if no hypothesis, set first coefficient to zero
    if (is.null(hypo)) hypo <- c(0, rep(1, length=(np1-1)))
    if (length(hypo) != np1) stop ("Hypothesis must be same length as no. of coefficients")
    m1 <- "Hypothesis must be a logical vector of TRUE/FALSE (or 0/1) statements"
    if (!all(c(0, 1) %in% hypo | !all(c(TRUE, FALSE) %in% hypo))) stop (m1)
### partial log-likelihood from original model
    lr1 <- x$loglik[2]
### taken from code for lm()
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
    if(!all) {
        r1 <- .getLR(DT=dt1, lr=lr1, hypo=hypo,
                     mFrame=mf, ties=ties, x=x)
### r1[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
        return(r1)
    }
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
### swap 1s with 0s to show which coefficients have been preserved
### i.e. keep those which are thought to be ==0
### res1 - holds results
    colnames(c2) <- n1
    lrdt1 <- data.table(c2)
### initialize with arbitrary values giving column type
    lrdt1[, c("chiSq", "df", "pVal") := list(0.1, 1, 0.1)]
###
    for (i in seq(nrow(c2))){
        hypo <- c2[i, ]
        lrdt1[i, ] <- .getLR(DT=dt1, lr=lr1, hypo=hypo,
                             mFrame=mf, ties=ties, x=x)
    }
    lrdt1 <- lrdt1[order(lrdt1$pVal), ]
### lrdt1[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
    return(lrdt1)
}
###
.getLR <- function(DT, lr, hypo, mFrame, ties, x){
### for R CMD check
    chiSq <- n1 <- NULL
### copy to modify (drop terms later)
    dt2 <- copy(DT)
    n1 <- names(x$coefficients)
    n2 <- n1[as.logical(hypo)]
    cn1 <- seq_along(colnames(dt2))
    i1 <- as.integer(cn1[!colnames(dt2) %in% n2])
    set(dt2, j=i1, value=NULL)
    stopifnot(attr(model.response(mFrame), "type")=="right")
    y1 <- data.table(unclass(model.response(mFrame, "numeric")))
    dt2[, c("t", "e") := y1]
### refit with only those coefficients
    lr2 <- survival::coxph(Surv(t, e) ~ .,
                           data=dt2,
                           ties=ties
                           )$loglik[2]
    r1 <- data.table(t(hypo),
                     2 * (lr - lr2),
                     sum(!hypo)
                     )
    setnames(r1, c(n1, "chiSq", "df"))
    r1[, "pVal" := (1-stats::pchisq(chiSq, df))]
    return(r1)
}
##' @rdname local
##' @aliases locWald
##' @export 
##'
locWald <- function(x, ...){
    UseMethod("locWald")
}
##' @rdname local
##' @aliases locWald.coxph
##' @method locWald coxph
##' @export
##' 
##' @examples
##' ###
##' data(larynx, package="KMsurv")
##' c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx, method="breslow")
##' locWald(c1, all=TRUE)
##' locWald(c1, hypo=c(0, 0, 0, 1))
locWald.coxph <- function(x, ...,
                          all=FALSE,
                          hypo=NULL
                          ){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### get names of coefficients
    n1 <- names(x$coefficients)
### no. predictor terms
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
### if no hypothesis, set first coefficient to zero
    if (is.null(hypo)) hypo <- c(0, rep(1, length=(np1-1)))
    if (length(hypo) != np1) stop ("Hypothesis must be same length as no. of coefficients")
    m1 <- "Hypothesis must be a logical vector of TRUE/FALSE (or 0/1) statements"
    if (!all(c(0, 1) %in% hypo | !all(c(TRUE, FALSE) %in% hypo))) stop (m1)
### one hypothesis only
    if(!all){
        r1 <- .getWald(x=x, hypo=hypo)
### r1[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
        return(r1)
    }
###
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
    colnames(c2) <- n1
### swap 1s with 0s to show which coefficients have been preserved
### i.e. keep those which are thought to be ==0
### w1 - holds results
    w1 <- data.table(1-c2)
    w1[, c("chiSq", "df", "pVal") := 0.1]
###
    for (i in seq(nrow(c2))){
        hypo <- c2[i, ]
        w1[i, ] <- .getWald(x=x, hypo=hypo)
    }
    w1 <- w1[order(w1$pVal), ]
### w1[, c("chiSq", "pVal") := list(formatC(chiSq), formatC(pVal))]
    return(w1)
}
###
.getWald <- function(x, hypo){
### for R CMD check
    chiSq <- n1 <- NULL
### coefficients to keep
### i.e. keep those which are thought to be ==0
    n1 <- names(x$coefficients)
    pos1 <- which(hypo==1)
    chi1 <- x$coefficients[pos1] %*% solve(x$var[pos1, pos1]) %*%  x$coefficients[pos1]
    df1 <- sum(hypo)
    r1 <- data.table(t(hypo),
                     chi1,
                     df1)
    setnames(r1, c(n1, "chiSq", "df"))
    r1[, "pVal" := (1-stats::pchisq(chiSq, df))]
    return(r1)
}
