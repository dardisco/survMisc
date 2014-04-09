###
##' @name genDfSurv
##' @export
##' @title Generate data frame for survival analysis
##' @description
##' Generates a data frame with a binary outcome (events), times and optional
##' predictor variables. An optional coxph model is fitted to the data.
##' Model is fitted with \code{survival::coxph()}
##'
##' @param b \dfn{binomial predictors}, the number of predictors which are
##' binary, i.e. limited to 0 or 1
##' @param f \dfn{factors}, the number of predictors which are factors
##' @param c \dfn{continuous predictors}, the number of predictors which are
##' continuous
##' @param n number of observations in the data frame
##' @param nf the no. of levels in a factor
##' @param pb \dfn{probability for binomnial predictors}
##' the probability of binomial predictors being =1.
##' E.g. if \code{pb=0.3}, 30\% will be 1s, 70\% will be 0s
##' @param rc \dfn{ratio for continuous variables} the ratio of levels of
##' continuous variables to the total number of observations \dfn{n} e.g. if
##' \code{rc=0.8} and \code{n=100}, it will be in the range 1-80
##' @param pe \dfn{probability of event}
##' the probability of events (typically death/failure) occurring,
##' i.e. are =1.
##' E.g. if \code{pe=0.5}, 50\% will be 1s, 50\% will be 0s
##' @param t0 Lowest (starting) time
##' @param tmax Highest (final) time
##' @param asFactor if \code{TRUE}, predictors given as factors
##' will be converted to factors in the data frame before the model is fit
##' @param timelim function will timeout after \code{timelim} secs
##' @param model If \code{TRUE} will return fitted model
##' @return A list with the following values:
##'  \item{df}{data frame with predictors
##' (labelled \eqn{x1,x2, ..., xn}) and outcome
##'   (y), with \emph{n} rows (observations)}
##'  \item{model}{if \code{model = TRUE} model fit with
##'   \code{ survival::coxph()} }
##' @note Using \code{asFactor=TRUE} with factors which have a large number of
##' levels (e.g. \code{nf >30}) on large datasets (e.g. \eqn{n >1000}) can cause
##' fitting to be slow.
##'
##' @keywords datagen
##' @examples
##' set.seed(1)
##' genDfSurv()
##' genDfSurv(b=0, c=2, n=100, pe=0.7)
##' genDfSurv(b=1, c=0, n=1000)
##' genDfSurv(f=1, nf=4, b=1, c=0, asFactor=TRUE)
genDfSurv <- function(f=0, b=2, c=1, n=100,
                      pb=0.5, nf=3, rc=0.8, pe=0.5,
                      t0=1, tmax=100,
                      asFactor=FALSE, model=TRUE, timelim=5) {
    if ( (pb|rc|pe) > 1 ) stop("Ratios should be <1")
    if (nf <= 2) stop("Factors should have at least 3 levels")
    bcf <- c(b, c, f)
    if (sum(bcf) <= 0) stop("Need at least one predictor")
    if (sum(bcf) >= n) stop("Need n to be larger for this no. predictors")
### prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
### define frame to hold values
    df1 <- as.data.frame(matrix(NA, ncol=(sum(bcf)+2), nrow=n))
    xnam <- paste("x", 1:(ncol(df1)-2), sep="")
    colnames(df1) <- c(xnam, "t1", "e")
### repeat until no dupicated columns
    repeat{
        if(!f==0){
            df1[, 1:f] <- sample(x=seq(1, nf),
                                 size=f*n,
                                 replace=TRUE)
        }
        if(!b==0){
            df1[, (f+1):(f+b)] <- sample(x=c(0L, 1L),
                                         size=b*n,
                                         replace=TRUE,
                                         prob=c(pb, 1-pb)
                                         )
        }
        if(!c==0){
            df1[, (b+f+1):(b+f+c)] <- sample(x=seq(1, (rc*n)),
                                             size=c*n,
                                             replace=TRUE)
        }
        df1[, ncol(df1)] <- sample(x=c(0L, 1L),
                                   size=n,
                                   replace=TRUE,
                                   prob=c(pe, 1-pe)
                                   )
        df1[, (ncol(df1)-1)] <- seq(from=t0, to=tmax, length.out=n)
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(df1))){break}
    }
### convert to factor
    if (!f==0 && asFactor==TRUE){
        df1[, 1:f] <- do.call(as.factor, list(df1[, 1:f] ))
        }
### generate model
    if(!model) return(list(df=df1))
###
    s1 <- survival::Surv(df1$t1, df1$e)
    f1 <- as.formula(paste("s1 ~ ", paste(xnam, collapse= "+")  ))
    c1 <- survival::coxph(f1, data=df1)
###
    res <- list(df=df1,
                model=c1)
    return(res)
}
