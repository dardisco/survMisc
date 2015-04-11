##' @name genSurv
##' @rdname genSurv
##' @aliases genSurvDf
##' @aliases genSurvDt
##' @title Generate survival data
##' 
##' @param b \dfn{binomial predictors}, the number of predictors which are
##' binary, i.e. limited to \eqn{0} or \eqn{1}
##' @param f \dfn{factors}, the number of predictors which are factors
##' @param c \dfn{continuous predictors}, the number of predictors which are
##' continuous
##' @param n number of observations (rows) in the data
##' @param nlf the number of levels in a factor
##' @param pb \dfn{probability for binomnial predictors}:
##' the probability of binomial predictors being \eqn{=1}
##' e.g. if \code{pb=0.3}, \eqn{30\%} will be \eqn{1}s, \eqn{70\%} will be \eqn{0}s
##' @param rc \dfn{ratio for continuous variables}: the ratio of levels of
##' continuous variables to the total number of observations \dfn{n} e.g. if
##' \code{rc=0.8} and \code{n=100}, it will be in the range \eqn{1-80}
##' @param pe \dfn{probability of event}
##' the probability of events (typically death/failure) occurring,
##' i.e. \eqn{P(e=1)}.
##' e.g. if \code{pe=0.6}, \eqn{60\%} will be \eqn{1}s, \eqn{40\%} will be \eqn{0}s
##' @param t0 Lowest (starting) time
##' @param tMax Highest (final) time
##' @param asFactor if \code{asFactor=TRUE} (the default),
##' predictors given as \code{factor}s
##' will be converted to \code{factor}s in the data frame before the model
##' is fit
##' @param timelim function will timeout after \code{timelim} secs.
##' This is present to prevent duplication of rows.
##' @param model If \code{model=TRUE} will also
##' return a model fitted with \code{survival::coxph}.
##' @return If \code{model=FALSE} (the default) a \code{data.frame} or \code{data.table} as above.
##' \cr
##' If \code{model=TRUE}: a list with the following values:
##'  \item{df or dt}{A \code{data.frame} (for \code{genSurvDf}) or \code{data.table}
##' (for \code{genSurvDt}).
##' Predictors are labelled \eqn{x1, x2, ..., xn}.
##' Outcome is \eqn{t1} (time) and \eqn{e} event (e.g. death).
##' Rows represent to \eqn{n} observations}
##'  \item{model}{A model fit with \code{survival::coxph}}
##' 
##' @note \code{genSurvDt} is faster and more efficient for larger datasets.
##' \cr \cr
##' Using \code{asFactor=TRUE} with factors which have a large number of
##' levels (e.g. \code{nlf >30}) on large
##' datasets (e.g. \eqn{n >1000}) can cause
##' fitting to be slow.
##' 
##' @keywords datagen
##' @keywords survival
##' 
##' @examples
##' set.seed(1)
##' genSurvDf(model=TRUE)
##' genSurvDf(b=0, c=2, n=100, pe=0.7)
##' genSurvDf(b=1, c=0, n=1000)
##' genSurvDf(f=1, nlf=4, b=1, c=0, asFactor=FALSE)
##' set.seed(1)
##' genSurvDt()
##' genSurvDt(b=0, f=0, c=1, n=20L, pe=0.7)
NULL
##'
##' @rdname genSurv
##' @export
##' 
genSurvDf <- function(b=2L, f=2L, c=1L, n=100L,
                      pb=0.5, nlf=3L, rc=0.8, pe=0.5,
                      t0=1L, tMax=100L,
                      asFactor=TRUE,
                      model=FALSE,
                      timelim=5) {
    if ( (pb|rc|pe) > 1 ) stop("Ratios should be <1")
    if (nlf <= 2) stop("Factors should have at least 3 levels")
    sbcf1 <- sum(c(b, f, c))
    if (sbcf1 <= 0) stop("Need at least one predictor")
    if (sbcf1 >= n) stop("Need n to be larger for this no. predictors")
### prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
### define frame to hold values
    df1 <- as.data.frame(matrix(NA, ncol=(sbcf1+2L), nrow=n))
    xnam1 <- paste("x", 1:(ncol(df1)-2), sep="")
    colnames(df1) <- c(xnam1, "t1", "e")
### repeat until no dupicated columns
    repeat{
        if(!b==0){
            df1[, 1:b] <- sample(x=c(0L, 1L),
                                         size=b*n,
                                         replace=TRUE,
                                         prob=c(pb, 1-pb)
                                         )
        }
        if(!f==0){
            if (asFactor){
            df1[, (b+1L):(b+f)] <- as.factor(sample(x=seq(1, nlf),
                                           size=f*n,
                                           replace=TRUE)
                                    )
            } else {
                df1[, (b+1L):(b+f)] <- sample(x=seq(1, nlf),
                                              size=f*n,
                                              replace=TRUE)
            }
        }
        if(!c==0){
            df1[, (b+f+1L):(b+f+c)] <- sample(x=seq(1, (rc*n)),
                                             size=c*n,
                                             replace=TRUE)
        }
        df1[, ncol(df1)] <- sample(x=c(0L, 1L),
                                   size=n,
                                   replace=TRUE,
                                   prob=c(pe, 1-pe)
                                   )
        df1[, (ncol(df1)-1)] <- seq.int(from=t0, to=tMax, length.out=n)
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(df1))){break}
    }
### generate model
    if(!model) return(df1)
###
    s1 <- survival::Surv(df1$t1, df1$e)
    f1 <- as.formula(paste("s1 ~ ", paste(xnam1, collapse= "+")  ))
    c1 <- survival::coxph(f1, data=df1)
###
    res <- list(df=df1,
                model=c1)
    return(res)
}
##'
##' @rdname genSurv
##' @export
##' 
genSurvDt <- function(b=2L, f=2L, c=1L, n=100L,
                      pb=0.5, nlf=3L, rc=0.8, pe=0.5,
                      t0=1L, tMax=100L,
                      asFactor=TRUE, model=TRUE, timelim=5) {
    if ( (pb|rc|pe) > 1 ) stop("Ratios should be <1")
    if (nlf <= 2) stop("Factors should have at least 3 levels")
    sbcf1 <- sum(c(b, f, c))
    if (sbcf1 <= 0) stop("Need at least one predictor")
    if (sbcf1 >= n) stop("Need n to be larger for this no. predictors")
### prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
### define frame to hold values
    dt1 <- data.table::data.table(matrix(0, ncol=(sbcf1+2L), nrow=n))
    xnam1 <- paste("x", 1:(ncol(dt1)-2), sep="")
    data.table::setnames(dt1, c(xnam1, "t1", "e"))
### repeat until no dupicated columns
### repeat until no dupicated columns
    repeat{
        if(!b==0){
            for (i in 1L:b){
                data.table::set(dt1, j=i, value=sample(x=c(0L, 1L),
                                          size=n,
                                          replace=TRUE,
                                          prob=c(pb, 1-pb)
                                          )
                                )
            }
        }
        if(!f==0){
            for (i in (b+1L) : (b+f)){
                if (asFactor){
                    data.table::set(dt1, j=i, value=as.factor(sample(x=seq.int(1L, nlf),
                                              size=n,
                                              replace=TRUE))
                                    )
                } else {
                    data.table::set(dt1, j=i, value=sample(x=seq.int(1L, nlf),
                                              size=n,
                                              replace=TRUE)
                                    )
                }
            }
        }
        if(!c==0){
            for (i in (b+f+1L) : (b+f+c)){
                data.table::set(dt1, j=i, value=sample(x=seq.int(1, (rc*n)),
                                          size=n,
                                          replace=TRUE)
                                )
            }
        }
        data.table::set(dt1, j=as.integer(sbcf1+1L), value=  seq(from=as.integer(t0),
                                         to=as.integer(tMax),
                                         length.out=n))
        data.table::set(dt1, j=as.integer(sbcf1+2L), value= sample(x=c(0L, 1L),
                                         size=n,
                                         replace=TRUE,
                                         prob=c(pe, 1-pe)
                                         ))
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(dt1))){break}
    }
### generate model
    if(!model) return(dt1)
###
    s1 <- survival::Surv(dt1$t1, dt1$e)
    f1 <- as.formula(paste("s1 ~ ", paste(xnam1, collapse= "+")  ))
    c1 <- survival::coxph(f1, data=dt1)
###
    res <- list(dt=dt1,
                model=c1)
    return(res)
}
