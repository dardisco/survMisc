##' @name profLik
##' @export
##' @title Profile likelihood for coefficients in Coxph model
##' @param x A \code{coxph} model
##' @param CI Confidence Interval
##' @param interval No. of points over which to evaluate coefficient.
##' @param mult Multiplier. Coefficent will be multiplied by lower and upper
##' value and evaluated across this range
##' @param noPerPage Number of plots per page (for initial plots).
##' Will be used as guidance and optimised for ease of display
##' @param width Width of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @param height Height of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @details
##' Plots of range of values for coefficient in model with log-likelihoods
##' for the model with the coefficient fixed at these values.
##' \cr \cr
##' For each coefficient a range of possible values is chosen, given by \eqn{\hat{B}*mult_{lower} -
##' \hat{B}*mult_{upper}}{Bhat*mult[lower] - Bhat*mult[upper]}.
##' A series of model are fit (given by \code{interval}). The coefficient is included in the model as a
##' \emph{fixed} term and the partial log-likelihood for the model is calculated.
##' \cr \cr
##' A curve is plotted which gives the partial log-likelihood for each of these candidate values.
##' An appropriate confidence interval (CI) is given by subtracting 1/2 the value of the approapriate quantile
##' of a chi-squared distribution with 1 degree of freedom.
##' \cr \cr
##' Two circles are also plotted giving the 95% CI for the Wald statistic.
##' @return One plot for each coefficient in the model
##' @examples
##' c1 <- coxph(formula = Surv(time, status == 2) ~ age + edema + log(bili) +
##' log(albumin) + log(protime), data = pbc)
##' profLik(c1)
##' @references Example is from:
##' Therneau T, Grambsch P 2000.
##' \emph{Modeling Survival Data}, 1st edition.
##' New York: Springer.
##' Section 3.4.1, pg 57.
##'
profLik <- function(x, CI=0.95, interval=50, mult=c(0.1, 2),
                    noPerPage=2, width=1500, height=800){
    if(!inherits(x, "coxph")) stop ("Only applies to objects of class coxph")
    coef1 <- stats::coef(x)
    l1 <- length(coef1)
###----------------------------------------
### set up window to plot
###
### find best balance for noPerPage
    if (noPerPage==1) {
        nrow1 <- ncol1 <- 1
        } else {
            balance <- function(x) (x^2 + noPerPage )/ x
            nrow1 <- round(stats::optimize(balance,
                                           interval=seq(1:noPerPage))$minimum,
                           digits=0)
            ncol1 <- round(noPerPage/nrow1, digits=0)
        }
### open plot window + save parameters before altering
    windows(record=TRUE, width=width, height=height)
    p <- graphics::par
### oma=outer margins, mar=margins, bottom,left,top,right
    graphics::par( mfrow=c(nrow1,ncol1), oma=c(0,0,4,0))
### plot title
    f1 <- as.character(stats::formula(x))
###  substitute >=2 whitespaces with 1
    f1 <- sub("  +", " ", x=f1)
    maintext <- function(){
        tex1 <- paste0("Partial likelihood profiles and ",
                       100*CI,
                       "% CI cutoff for model ",
                       f1,
                       " \n Circles show ",
                       100*CI,
                       "% CI limits for Wald interval",
                       sep="")
        graphics::mtext(tex1, line=0.3, outer=TRUE)
        }
###
###----------------------------------------
###  plots
###
### get names of the coefficients from model.frame
### note excluding Surv
    n1 <- names(model.frame(x))[!grepl( "Surv", names(model.frame(x)) )]
### allocate memory for
    llik <- double(length=interval)
###
    for (i in 1:l1){
### lower + upper limits
        low1 <- mult[1]*coef1[i]
        up1 <- mult[2]*coef1[i]
### range for coefficient
        beta1 <- seq(from=low1, to=up1, length.out=interval)
        for (j in 1:interval){
### right hand side of formula without coefficient
            rhs1 <- paste0(n1[-i], collapse="+")
### offset = includes coefficient as fixed covariate
            off1 <- beta1[j]
            off2 <- paste0("+offset(",off1,"*",n1[i],")")
### new RHS for formula
            rhs2 <- paste0(rhs1, off2)
            f2 <- stats::as.formula(paste0(".~", rhs2))
### refit model and find model loglik with this value (beta) of coefficient
            c2 <- stats::update(x, formula=f2)
            llik[j] <- c2$loglik[2]
            }
        graphics::plot(beta1, llik, type="l",
             xlab="Values for coefficient",
             ylab="Model partial likelihood",
             main=n1[i]
             )
### range for confidence interval is chi-square on with 1 df
        rCI <- stats::qchisq(CI,1)
### confidence interval (calcuate lower only)
        ci1 <- x$loglik[2] -rCI/2
        graphics::abline(h=ci1, lty=2)
        sd1 <- sqrt(x$var[i,i])
### range for confidence interval of Wald is normal
### if CI is 95% then need convert to 97.5%
        CI2 <- (1-CI)/2
        rCI <- stats::qnorm(1-CI2)
        graphics::points(coef1[i] +c(-rCI, rCI)*sd1, c(ci1, ci1),
                         pch=1, cex=3)
        maintext()
    }
    par <- p
}
