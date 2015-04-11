##' @name sig
##' @title Significiance tests of coefficients in a \code{coxph} model
##'
##' @keywords htest
##' 
##' @rdname sig
##' @export
sig <- function(x, ...){
    UseMethod("sig")
}
##'
##' @rdname sig
##' @aliases sig.coxph
##' @aliases log-rank test
##' @method sig coxph
##' @export
##'
##' @param x A model of class \code{coxph}
##' @param ... Additional arguments (not implemented)
##' 
##' @return A \code{data.frame} with one row for each coefficient in the
##' original model. There are three columns, one for each of the tests:
##' \item{Wald}{the statistic is:
##'             \deqn{\frac{\hat{B}}{\hat{SE}}}{
##'                   B/SE}
##'             where \eqn{\hat{B}}{B} is the estimate of the coefficient
##'             and \eqn{\hat{SE}}{SE} is its standard error.}
##' \item{plr}{\bold{P}artial \bold{l}ikelihood \bold{r}atio test.
##'            \cr
##'           The statistic is the difference in the
##'           likelihood ratio of the original model and that with the coefficient
##'           omitted.}
##' \item{lrt}{Aka the \bold{score} test.
##'              \cr
##'              The Null hypothesis is that
##'              \eqn{\hat{B}=0}{B=0}.
##'              \cr
##'              The statistic is cacluated by refitting the model with the coefficient
##'              omitted, to generate initial values for the other \eqn{\hat{B}}{B}s.
##'              \cr
##'              It is then fitted again with all
##'              covariates, using these values and setting \eqn{\hat{B}=0}{B=0}.
##' }
##' All statistics are distributed as \eqn{\chi}{chi}-square, with degrees of freedom
##' \eqn{=} no. of coefficients \eqn{-1}.
##'
##' 
sig.coxph <- function(x, ...){
    if(!inherits(x, "coxph")) stop
    ("Only applies to objects of class coxph")
    l1 <- length(coefficients(x))
    if (l1==0) stop
    ("No coefficients; this is an intercept-only model")
### hold results
    res1 <- data.frame(matrix(0.01, nrow=l1, ncol=3))
### std. errors
    se1 <- sqrt(diag(x$var))
### p value for Wald tests
    res1[ ,1] <- 1 - pchisq((coef(x)/se1) ^ 2, 1)
### likelihood ratio test statistic
    LR1 <- -2 * (x$loglik[1] - x$loglik[2])
### get names of the coefficients from model.frame
### note excluding Surv
    n1 <- colnames(model.matrix(x))
### if only one coefficient then will be vs intercent-only model
    if (l1==1){
        res1[1, 2] <- 1 - pchisq(LR1, 1)
        res1[1, 3] <- 1 - pchisq(x$score, 1)
        rownames(res1) <- n1
        colnames(res1) <- c("Wald", "LR", "score")
        return(res1)
    }
### find degrees of freedom ( taken from survival:::print.coxph() )
    findDf <- function(x){
        if (is.null(x$df)) {
            df <- sum(!is.na(x$coefficients))
        } else {df <- round(sum(x$df), 2)}
    }
    degf1 <- findDf(x)
### log-likelihood for original model
    pLL <- 1 - stats::pchisq(LR1, degf1)
    for (i in 1:length(n1)){
### ### refit with coefficient omitted
        y1 <- model.response(model.frame(x))
        x1 <- model.matrix(x)[, -i]
        c2 <- coxph(y1 ~ x1)
        degf2 <- findDf(c2)
        LR2 <- -2 * (c2$loglik[1] - c2$loglik[2])
        LRdiff <- LR1 - LR2
        dfDiff <- degf1 - degf2
        pLL <- 1 - pchisq(LRdiff, dfDiff)
        res1[i, 2] <- pLL
        if(!x$n==c2$n) warning
        ("Need same no. observations to compare models; check for missing data ")
        inits1 <- coef(x)
        inits1[i] <- 0
### ### refit with new initial coefficients
### ### no iterations 
        c3 <- update(x, init=c(inits1), iter=0)
        res1[i, 3] <- pchisq(c3$score, dfDiff)
    }
    rownames(res1) <- n1
    colnames(res1) <- c("Wald", "plr", "lrt")
    return(res1)
}
