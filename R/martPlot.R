##' @name martPlot
##' @include multPlot.R
##' @export martPlot
##' @title Martingale residuals from null model vs coefficients in Coxph model
##' @param x A \code{coxph} model
##' @param noPerPage Number of plots per page
##' Will be used as guidance and optimised for ease of display
##' @param width Width of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @param height Height of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @details
##' The Cox model assumes that the hazard function satisfies:
##' \deqn{\lambda_{i}(t) = \lambda_0(t) \exp(X_i\beta)}{
##'  h[i](t) = h[0](t) exp (X[i]Beta)}
##' That is, for a continuous variable, an unit increase in the variable produces the same
##' change in risk across the value of the variable. (E.g. an increase of age of 5 years leads to the same
##' change in hazard no matter what the increase is from or to).
##' \cr \cr
##' To verify this is the case, a null model is fitted
##' (i.e no coefficients, similar to intercept-only model in linear regression).
##' Martingale residuals are calcuated for this. Plots of these residuals against the values of each of the
##' predictors in the model are shown.
##' \cr \cr
##' If the correct model for covariate \eqn{j} is based on a smooth function \eqn{f()}, i.e.
##' \eqn{ \exp(f(X_j)\beta_j)}{exp(f(X[j])Beta[j])} then the following should hold:
##' \deqn{E(M_i|X_{ij}=x_j) \approx cf(x_j)}{
##'  E( M[i]| X[i,j]=x[j] ) = cf(x[j])}
##' Where \eqn{M} is the martingale residual and \eqn{c} depends on the amount of censoring and is
##' roughly independent of \eqn{x_j}{x[j]}.
##' \cr \cr
##' A \code{lowess} smoothed line is added to the plot. This should be approximately linear if the
##' assumption of proportional hazards is met.
##' @return One plot per coefficient, showing the values for the vaiable against the martingale
##' residuals from the null model.
##' @examples
##' c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
##' martPlot(c1)
##' @references Example is from:
##' Therneau T, Grambsch P 2000.
##' \emph{Modeling Survival Data}, 1st edition.
##' New York: Springer.
##' Section 5.1.2, pg 91.
##'
martPlot <- function(x, noPerPage=2, width=1500, height=800){
    if(!inherits(x, "coxph")) stop ("Only applies to objects of class coxph")
### set up for plots
    .multPlot(noPerPage=noPerPage, width=width, height=height)
### plot title
    f1 <- as.character(stats::formula(x))
###  substitute >=2 whitespaces with 1
    f1 <- sub("  +", " ", x=f1)
    maintext <- function(){
        tex1 <- paste0("Martingale Residuals for coefficients in model ",
                       f1,
                       "\n (should be linear if variable is linear in model)",
                       sep="")
        graphics::mtext(tex1, line=0.3, outer=TRUE)
        }
### refit as intercept only
    f2 <- stats::as.formula(".~1")
    c2 <- stats::update(x, formula=f2)
### get residuals (hidden method in survival)
    r1 <- survival:::residuals.coxph(c2, type="martingale")
### get names of the coefficients from model.frame
### note excluding Surv
    n1 <- names(model.frame(x)[ , -1, drop=FALSE])
    for (i in 1:length(n1)){
        x1 <- model.frame(x)[ , -1, drop=FALSE][ ,i]
        plot(x1, r1,
             xlab=n1[i],
             ylab="Residual = excess events (deaths)")
### add smoothing line
### use iter=0 to prevent removal of potential outliers
        graphics::lines(stats::lowess(x1, r1, iter=0), lty=1)
        maintext()
    }
}
