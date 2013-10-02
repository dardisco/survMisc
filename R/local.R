##' @name local
##' @rdname local
##' @title Local tests for a model
##' @export locScore
##' @aliases locScore
locScore <- function(x, ...){
    UseMethod("locScore")
}
##' @rdname local
##' @aliases locScore.coxph
##' @method locScore coxph
##' @S3method locScore coxph
##'
##' @include renameFact.R
##' @include getSurvNames.R
##'
##' @param x A model of class \code{coxph}
##' @param ... Additional arguments
##' @param ties Method of handling ties when refitting model. Must be
##' one of \code{breslow} or \code{efron}
##' @param all Fit \emph{all} combinations of predictors
##' @param hypo Hypothesis to test. This should be a vector of
##' \eqn{0}s and \eqn{1}s of the same length as the no. of
##' coefficients in the model (with at least one 0 and one 1).
##' Zeros indicate coefficients to exclude, ones indicate coefficients
##' to keep.
##' \cr
##' @return For \code{locScore} a \code{list} with the following items:
##' \cr
##' \item{coefficients}{coefficients from refitted model(s)}
##' \item{test}{test for hypothesis}
##' \cr \cr
##' For \code{locLR} and \code{locWald}, a \code{data.frame}
##' showing the hypothesis being tested and the results of the test.
##'
##' @details The null hypothesis is that some of the coefficients in
##' the model are zero (\eqn{H_0 :  \hat{B}_i=0, \quad i \geq 1}{H0: Bi=0})
##' vs. the alternative that at least one of them is nonzero.
##' \cr \cr
##' All of these tests are distributed as chi-square
##' with degrees of freedom = number of excluded coefficients.
##' \cr \cr
##' For the \bold{score} test,
##' the model is fitted again with the coefficients of
##' interest excluded.
##' \cr
##' A value for the remaining coefficients is obtained.
##' Then the complete model is fit again using these new values as initial values for
##' those remaining coefficients and using zero as the initial value for the excluded coefficients.
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
##' locScore(c1, hypo=c(0,0,0,1))
##' @references Examples are from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 8.2, pp 264-6.
locScore.coxph <- function(x, ...,
                           all=FALSE,
                           hypo="",
                           ties="breslow"){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### predictor terms
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
    n1 <- names(x$coefficients)
### if no hypothesis, use alternating coefficients
    if (hypo[1]=="") hypo <- rep(c(0,1), length.out=np1)
    if (!length(hypo)==np1) stop ("Hypothesis must be same length as no. of coefficients")
    if (!all(c(0,1) %in% hypo)) stop ("Hypothesis must be 0 AND 1 only")
### get names for time and event
### from Surv object in formula
    y1 <- .getSurvNames(x)
### one hypothesis only
    if(!all){
### coefficients to keep
        x1 <- n1[as.logical(hypo)]
### rename factor coefficients
        x1 <- .renameFact(x1)
### coefficients == 0
        x2 <- n1[!hypo]
        x2 <- .renameFact(x2)
### refit with only those coefficients
        f1 <- stats::as.formula(
            paste("Surv(", y1[1], ", ", y1[2], ") ~ ",
                  paste(x1, collapse="+"), sep="") )
        data1 <- eval(x$call$data)
        coef1 <- survival::coxph(f1,
                                 data=data1,
                                 ties=ties
                                 )$coefficients
### if 2x coefficients for 1 term, keep only those == TRUE
### eg x1:factor(x2==a)TRUE (vs FALSE)
        if (any(grepl("TRUE",names(coef1)))){
            coef1 <- coef1[!grepl("FALSE",names(coef1))]
        }
### add to results
        coef1 <- round(coef1, 3)
        init1 <- hypo
        init1[which(hypo==1)] <- coef1
        init1 <- as.vector(init1)
### get original formula
### need deparse then as.formula to reset environment
### in which formula called
        f2 <- stats::as.formula(deparse(x$formula))
### refit with new initial values and no iterations
### (i.e. no convergence of estimates)
        cox1 <- survival::coxph(f2,
                                data=data1,
                                ties=ties,
                                init=init1,
                                iter.max=0)
### score vector
        sc1 <- colSums(survival::coxph.detail(cox1)$score)
### chi-square
        chiSq <- round(sc1 %*% cox1$var %*% sc1, 3)
        df <- sum(hypo==0)
        pVal <- round(1-stats::pchisq(chiSq, df), 5)
        res1 <- data.frame(chiSq, df, pVal)
        dim(hypo) <- c(1,np1)
### result
        names(init1) <- colnames(hypo) <- n1
        return(list(coefficients=init1,
                    score=cbind(hypo,res1)))
    }
###
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
### drop first and last rows (i.e. all present or all absent)
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
    colnames(c2) <- n1
### add additional columns to store values
    res1 <- transform(c2, chiSq=0, df=0, pVal=0)
    colnames(res1)[1:length(n1)] <- n1
    y1 <- .getSurvNames(x)
    for (i in 1:nrow(c2)){
### coefficients to keep
        x1 <- n1[as.logical(c2[i, ])]
### rename factor coefficients
        x1 <- .renameFact(x1)
### coefficients == 0
        x2 <- n1[!c2[i, ]]
        x2 <- .renameFact(x2)
### refit with only those coefficients
        f1 <- stats::as.formula(
            paste("Surv(", y1[1], ", ", y1[2], ") ~ ",
                  paste(x1, collapse="+"), sep="") )
        data1 <- eval(x$call$data)
        coef1 <- survival::coxph(f1,
                                 data=data1,
                                 ties=ties
                                 )$coefficients
### if 2x coefficients for 1 term, keep only those == TRUE
### eg x1:factor(x2==a)TRUE (vs FALSE)
        if (any(grepl("TRUE",names(coef1)))){
            coef1 <- coef1[!grepl("FALSE",names(coef1))]
        }
### add to results
        c2[i, which(c2[i, ]==1)] <- round(coef1, 3)
### get original formula
        f2 <- stats::as.formula(deparse(x$formula))
        cox1 <- survival::coxph(f2,
                                data=data1,
                                ties=ties,
                                init=c2[i, ],
                                iter.max=0)
### score vector
        sc1 <- colSums(survival::coxph.detail(cox1)$score)
### chi-square
        res1$chiSq[i] <- round(sc1 %*% cox1$var %*% sc1, 3)
        df1 <- sum(c2[i, ]==0)
        res1$df[i] <- df1
        res1$pVal[i] <- round(1-stats::pchisq(res1$chiSq[i], df1), 5)
    }
###
    res1 <- res1[order(res1$pVal), ]
    c2 <- c2[order(res1$pVal), ]
    return(list(
        coefficients=c2,
        score=res1
        ))
}
##' @rdname local
##' @aliases locLR
##' @export locLR
##'
locLR <- function(x, ...){
    UseMethod("locLR")
}
##' @rdname local
##' @aliases locLR.coxph
##' @method locLR coxph
##' @S3method locLR coxph
##' @examples
##' ###
##' data(larynx, package="KMsurv")
##' c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx, method="breslow")
##' locLR(c1, all=TRUE)
##' locLR(c1, hypo=c(0,0,0,1))
##'
locLR.coxph <- function(x, ...,
                        all=FALSE, hypo="",
                        ties="breslow"){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### no. predictor terms
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
    n1 <- names(x$coefficients)
### if no hypothesis, use alternating coefficients
    if (hypo[1]=="") hypo <- rep(c(0,1), length.out=np1)
    if (!length(hypo)==np1) stop ("Hypothesis must be same length as no. of coefficients")
    if (!all(c(0,1) %in% hypo)) stop ("Hypothesis must be 0 AND 1 only")
### partial log-likelihood from original model
    lr1 <- x$loglik[2]
### get names for time and event
### from Surv object in formula
    y1 <- .getSurvNames(x)
### one hypothesis only
    if(!all){
### coefficients to keep
        x1 <- n1[as.logical(hypo)]
### rename factor coefficients
        x1 <- .renameFact(x1)
### refit with only those coefficients
        f1 <- stats::as.formula(
            paste("Surv(", y1[1], ", ", y1[2], ") ~ ",
                  paste(x1, collapse="+"), sep="") )
        lr2 <- survival::coxph(f1,
                               data=eval(x$call$data),
                               ties=ties
                               )$loglik[2]
        chiSq <- 2 * (lr1 - lr2)
        df <- sum(hypo==0)
        pVal <- round(1-stats::pchisq(chiSq, 3), digits=5)
        hypo <- t(as.matrix(hypo))
        colnames(hypo) <- n1
        res1 <- cbind(hypo, chiSq, df, pVal)
        colnames(res1)[(ncol(hypo)+1):ncol(res1)] <- c("chiSq", "df", "pVal")
        return(res1)
    }
###
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
### swap 1s with 0s to show which coefficients have been preserved
### i.e. keep those which are thought to be ==0
### res1 - holds results
    res1 <- transform(c2, chiSq=0, df=0, pVal=0)
###
    for (i in 1:nrow(res1)){
### coefficients to keep
        x1 <- n1[as.logical(c2[i, ])]
### rename factor coefficients
        x1 <- .renameFact(x1)
### refit with only those coefficients
        f1 <- stats::as.formula(
            paste("Surv(", y1[1], ", ", y1[2], ") ~ ",
                  paste(x1, collapse="+"), sep="") )
        lr2 <- survival::coxph(f1,
                               data=eval(x$call$data),
                               ties=ties
                               )$loglik[2]
        res1[i, "chiSq"] <- chiSq <- 2 * (lr1 - lr2)
        res1[i, "df"] <- sum(c2[i, ]==0)
        res1[i, "pVal"] <- round(1-stats::pchisq(res1[i, "chiSq"],
                                                 res1[i, "df"]),
                                 digits=5)
    }
    colnames(res1)[1:np1] <- n1
    return(res1)
}
##' @rdname local
##' @aliases locWald
##' @export locWald
##'
locWald <- function(x, ...){
    UseMethod("locWald")
}
##' @rdname local
##' @aliases locWald.coxph
##' @method locWald coxph
##' @S3method locWald coxph
##' @examples
##' ###
##' data(larynx, package="KMsurv")
##' c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx,
##' method="breslow")
##' locWald(c1, all=TRUE)
##' locWald(c1, hypo=c(0,0,0,1))
locWald.coxph <- function(x, ...,
                          all=FALSE, hypo=""){
    if (!class(x)=="coxph") stop("Only applies to objects of class coxph")
### no. predictor terms
    np1 <- length(x$coefficients)
    if (np1==1) stop ("Need >=2 coefficients")
    n1 <- names(x$coefficients)
### if no hypothesis, use alternating coefficients
    if (hypo[1]=="") hypo <- rep(c(0,1), length.out=np1)
    if (!length(hypo)==np1) stop ("Hypothesis must be same length as no. of coefficients")
    if (!all(c(0,1) %in% hypo)) stop ("Hypothesis must be 0 AND 1 only")
### one hypothesis only
    if(!all){
### coefficients to keep
### i.e. keep those which are thought to be ==0
        pos1 <- which(hypo==0)
        chiSq <- round(
            x$coefficients[pos1] %*%
            solve(x$var[pos1, pos1]) %*%
            x$coefficients[pos1],
            digits=3)
        df <- sum(hypo==1)
        pVal <- round(1-stats::pchisq(chiSq, df), digits=5)
        hypo <- t(as.matrix(hypo[,drop=FALSE]))
        colnames(hypo) <- n1
        res1 <- cbind(hypo,chiSq,df,pVal)
        colnames(res1)[(ncol(hypo)+1):ncol(res1)] <- c("chiSq", "df", "pVal")
        return(res1)
        }
###
### otherwise
### all combinations of dropped coefficients
    c2 <- combinat::hcube(rep(2, np1))-1
    c2 <- c2[-c(1, nrow(c2)), , drop=FALSE]
### swap 1s with 0s to show which coefficients have been preserved
### i.e. keep those which are thought to be ==0
### res1 - holds results
    res1 <- transform(1-c2, chiSq=0, df=0, pVal=0)
    colnames(res1)[1:length(n1)] <- n1
    for (i in 1:nrow(c2)){
        pos1 <- which(c2[i, ]==1)
        res1$chiSq[i] <- round(
            x$coefficients[pos1] %*%
            solve(x$var[pos1, pos1]) %*%
            x$coefficients[pos1],
            digits=3)
        df1 <- sum(c2[i,]==1)
        res1$df[i] <- df1
        res1$pVal[i] <- round(1-stats::pchisq(res1$chiSq[i],df1), digits=5)
        }
    res1 <- res1[order(res1$pVal), ]
    return(res1)
}
###
