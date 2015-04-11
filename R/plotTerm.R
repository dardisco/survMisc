##' @name plotTerm
##' @include gamTerms.R
##' @export
##' @title Plot individual terms of a Generalized Additive \code{gam}
##' or Cox Proportional Hazards \code{coxph} Model
##' 
##' @param x The result of a Generalized Additive Model (\code{gam})
##' or a Cox proportional hazards model (\code{coxph})
##' @param term The term to be plotted.
##' An integer, based on the position of the term on the right-hand side
##' of the model formula e.g. the first term is \code{1}. 
##' @param se If \code{se=TRUE} (the default), also plot confidence intervals
##' based on the standard errors
##' @param p P-value used to plot the confidence intervals
##' from standard errors. Based on the normal distribution
##' @param rug Add rug (1-dimensional plot) to x-axis.
##' See \code{?graphics::rug}
##' @param const Value for constant term. If \code{const=TRUE},
##' add the overall mean for the decomposition as returned by
##' \code{\link{gamTerms}}
##' @param col Color of line(s) on plot. If \code{se=TRUE}, use a
##' \code{vector} of 3 colors: the first is the main line,
##' the second is lower the CI, the third is the upper CI.
##' @param data A \code{data.frame} in which to evaluate the model.
##' If missing, \code{eval(x$call$data)} is used.
##' @param xlab Label for x-axis.
##' @param ... Additional arguments; passed to \code{graphics::matplot}
##' (with standard errors) or \code{graphics::plot} (without).
##' 
##' @return A plot (base graphics) of the term in question.
##' If \code{se=TRUE},
##' this is done using \code{graphics::matplot} otherwise
##' \code{graphics::plot.default} is used.
##' 
##' @author Terry Therneau. Updated from S-plus by Chris Dardis
##' @examples
##' fit1 <- coxph(Surv(time, status) ~ sex + pspline(age), data=lung)
##' plotTerm(fit1, term=2, rug=FALSE, ylab="Log-hazard",
##'          col=c("blue", "red", "red"))
##' ### smoothing splines
##' data(mpip)
##' c1 <- coxph(Surv(futime, status) ~
##'             pspline(lved) + factor(nyha) + rales + pspline(ef),
##'             data=mpip)
##' plotTerm(c1, 4, ylab="Log-hazard", xlab="Ejection fraction (%)",
##'          main="Log-hazard by ejection fraction \n Line fitted by penalized smoothing splines")

##' @seealso \code{\link{gamTerms}}, \code{\link{mpip}}
##' 
plotTerm <- function(x,
                     term=1,
                     se=TRUE,
                     p=0.95,
                     rug=TRUE,
                     const=0,
                     col=1,
                     xlab=NULL,
                     data,
                     ...) {
    stopifnot(inherits(x, "coxph") | inherits(x, "gam"))
    if(missing(data)) data <- eval(x$call$data)
    g1 <- gamTerms(x, data=data)
###----------------------------------------

    if(missing(xlab)) xlab <- names(g1)[term+1]
    xy1 <- g1[[term+1]]
    if (!is.numeric(const)) {
	if (const==TRUE) const <- g1[[1]]
    } else {
        const <- 0
    }
###
    if (se){
	ci <- -1 * qnorm((1 - p) / 2)
	ymat1 <- cbind(xy1[, "y"],
                      xy1[, "y"] - ci * xy1[, "se"],
                      xy1[, "y"] + ci * xy1[, "se"]) + const
	graphics::matplot(xy1[, "x"], ymat1,
                          type='l', ...,
                          xlab=xlab, lty=c(1, 2, 2),
                          col=col)
	if (rug) graphics::rug(xy1[, "x"])
    } else {
        graphics::plot(xy1[, "x"], xy1[, "y"] + const,
             type='l', ...,
             xlab=xlab, ylab="")
    }
}
