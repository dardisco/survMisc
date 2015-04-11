##' @name dx
##' @rdname dx
##' 
##' @title Diagnostics for \code{coxph} models
##' @export
##' 
dx <- function(x, ...){
    UseMethod("dx")
}

##' @rdname dx
##' @aliases dx.coxph
##' @method dx coxph
##' @export
##'
##' @include asSurv.R
##' @include plotTerm.R
##' @include gamTerms.R
##' @include sf.R
##' 
##' @param x An object of class \code{coxph}.
##' @param ... Additional arguments. Can be passed to \code{graphics::plot}
##' or \code{graphics::matplot}.
##' \cr
##' See \code{?par} for details.
##' @param what Which plots to make. See \bold{Value} below.
##' @param toPdf Print plots to pdf. This is usually recommended as
##' each plot is created on a new device and 'R' can typically only
##' have 61 devices open simultaneously.
##' \itemize{
##'   \item If \code{toPdf=TRUE}, each plot is created on a new page.
##'   \item If \code{toPdf=FALSE}, each plot is created on a new screen device.
##'  }
##' @param file Filename to store plots. Default is \code{"dxPlots.pdf"}.
##' @param maxStrata Used for time vs. log-log survival plot.
##' \cr
##' If there are \eqn{>} \code{maxStrata} strata, no plot is shown for this.
##' \cr
##' Recommended is \eqn{\leq 5}{<=5} to prevent the plot from becoming visually cluttered.
##' @param defCont Definition of continuous variable.
##' \cr
##' Variables with more than \code{defCont} unique values will
##' be split into quantiles to facilitate graphs.
##' \cr
##' (This does \emph{not} apply to \code{factor} variables).
##' @param noQuantiles No. of quantiles into which to
##' split continuous variables
##' @param maxFact Maximum number of levels in a factor.
##' \cr
##' Used in plotting differences in log-hazard curves.
##' @param identify Identify outliers manually. Cannot be used with \code{toPdf=TRUE}.
##' @param degfP Degrees of freedom for smoothing spline in Poisson model.
##' @param degfRS Degrees of freedom for regression splines.
##' @param degfSS Degrees of freedom for smoothing splines.
##' If \code{degFSS=0}, the 'optimal' degrees of freedom is chosen
##' according to the \code{\link{AIC}}.
##' @param timeTrans Type of time transformation to use when refitting model with
##' time-transformed covariate.
##' \cr
##' See \code{?survival::cox.zph} for details.
##' @param ties  Method of handling ties when refitting model (for stratified plots).
##' \cr
##' Default is the same as the original model, \code{x}.
##' Usually one of \code{"breslow"} or \code{"efron"}.
##'
##' @return
##' Plots with base graphics.
##' \cr
##' If \code{what="ph"}:
##'   \item{\eqn{\pm}{+-}}{Time vs. \eqn{-\log -\log}{-log-log} survival.
##'                          \cr
##'                          If not too many strata to plot, as per argument \code{maxStrata}.}
##'   \item{*}{Quantile-quantile plot.
##'            \cr
##'            Unit exponential distribution
##'            vs. expected events (or Cox-Snell residuals)}
##'   \item{*}{Observed vs. expected hazard}
##'   \item{*}{Expected events vs. hazard based on sorted expected events}
##'   \item{*}{Time vs. hazard, per predictor.
##'            \cr
##'            Continuous variables are split into quantiles.}
##'   \item{*}{Time vs. difference in log hazards, per predictor.
##'            \cr
##'            Continuous variables are split into quantiles.}
##'   \item{*}{Reference hazard vs. hazard for predictor.
##'            \cr
##'            Continuous variables are split into quantiles.}
##' If \code{what="lin"} (only applies to continuous variables):
##'   \item{\eqn{\pm}{+-}}{Predictor vs. residuals from a Poisson model with smoothing spline.}
##'   \item{\eqn{\pm}{+-}}{Predictor vs. partial residual for predictor (with regression spline).
##'                        \cr
##'                        For predictors with \eqn{>} \code{degfRS}. }
##'   \item{\eqn{\pm}{+-}}{Predictor vs. partial residual for predictor (with smoothing spline).
##'                        \cr
##'                        For predictors with \eqn{>} \code{degfSS}.}
##'   \item{*}{Time vs. scaled Schoenfeld residuals, per predictor.}
##' If \code{what="inf"}:
##'   \item{*}{Observation vs. jacknife influence.}
##'   \item{*}{Observation vs. jacknife influence (scaled by standard error of coefficients).}
##'   \item{*}{Observation vs. leverage (=scaled score residual).}
##'   \item{*}{Martingale residuals vs. likelihood displacement residuals.}
##' 
##' If \code{what="lin"}, a \code{list} of \code{data.table}s is also returned to the console:
##'
##'   \item{Poisson}{Results from anova for a Poisson fit (via \code{gam}) with nonparametric effects.
##'                  \cr
##'                  The model is re-fit with smoothing splines for continuous variables.
##'                  \cr
##'                  Needs at least one predictor to have \eqn{>3} unique values.}
##'   \item{tt}{Results from the time-transformed fit, using \code{survival::cox.zph}.}
##'
##' @note
##' \bold{TESTS OF PROPORTIONAL HAZARDS}
##' \cr \cr
##' A simple graphical test of proportional hazards, applicable to time-fixed
##' variables with a small number of levels, is a plot of
##' time vs. \eqn{-\log (-\log [\hat{S}(t)])}{-log(-log[S(t)])}.
##' \cr
##' The Kaplan-Meier curves should be parallel as:
##' \deqn{\hat{S}_i(t) = \exp -H_i(t) = \exp[-(H_0(t) \exp[\hat{\beta} X_i(t)])]}{
##'        S[i](t) = exp(-H[i](t)) = exp(- H[0](t).exp( B.X[i](t) ) )}
##' where \eqn{H_0(t)}{H[0](t)} is the Breslow's estimator of the baseline hazard
##' (i.e. all co-variates are zero),
##' often represented as \eqn{\lambda_0(t)}{lambda[0](t)}. Thus
##' \deqn{-\log (-\log [\hat{S}(t)]) = -\log(H_0(t)) - \hat{\beta} X_i(t)}{
##'       -log(-log[S(t)]) = -log(H[0](t)) - B.X[i]}
##'
##' \subsection{A note on Cox-Snell residuals}{
##'   Given \eqn{n} observations, the residuals are:
##'   \deqn{ \hat{M}_i = Y_i - E(Y_i), \quad i=1,...,n }{
##'          M[i] = Y[i] - E(Y[i]) }
##'   where \eqn{Y_i}{Y[i]} are the observed events,
##'   \eqn{E(Y_i)}{E(Y[i])} are the expected events and 
##'   \eqn{\hat{M}_i}{M[i]} is the vector of residuals, known as \bold{martingale} residuals.
##'   \cr
##'   The expected events \eqn{E(Y_i)}{E(Y[i])} are generated for each observation as
##'   \deqn{E(Y_i) = H_0(t) \exp \sum \hat{\beta} X_i(t)}{
##'         E(Y[i]) = H[0](t) exp( SUM B.X[i](t) )}
##'   The equation for these residuals may be rewritten as:
##'   \deqn{ E(Y_i) = Y_i - \hat{M}_i, \quad i=1,...,n }{
##'          E(Y[i]) = Y[i] - M[i], for i=1 to n}
##'   Somewhat unintuitively, these predicted values \eqn{E(Y_i)}{E(Y[i])}, are also known as the
##'   \bold{Cox-Snell} residuals.
##' }
##' 
##' These Cox-Snell residuals are used to assess the fit of a
##' proportional hazards model.
##' More formally, they are are generated from the (non time-dependent)
##' covariates \eqn{X}, a matrix
##' with one row per observation (total \eqn{n}) and additional indicators
##' of time \eqn{t} and status \eqn{\delta}{D}.
##' The estimated coefficients are \eqn{ \hat{\beta} }{B},
##' where \eqn{\hat{\beta}}{B} is a vector of length \eqn{p} (the number of predictors).
##' \cr
##' The residuals are:
##' \deqn{r_i = H_0(t_i) \exp \sum_{k=1}^p \hat{\beta_k}X_{ik}(t_i),\quad i=1,...,n,\quad k=1,...,p}{
##'        r[i] = H[0](t[i]) exp ( SUM B[k].X[ik](t[i]) ), for j=1 to n, k=1 to p}
##' If the coefficients are close to their true values, then \eqn{r_i}{r[i]}
##' should follow a unit-exponential distribution, i.e. \eqn{H_0(t) \approx t}{H[0](t) = t}.
##' \cr \cr
##' Thus a \code{qqplot} of \eqn{r_i}{r[i]} against a unit-exponential distribution 
##' is given. This is of limited value, as the \emph{null} model (no coefficients) will be
##' closer to the true exponential distribution than that with \emph{any} coefficients.
##' \cr \cr
##' Another simple graphical test is a plot of observed vs. expected values
##' for the cumulative hazard.
##' \cr
##' The expected values of the hazard are generated using the expected events
##' (using \code{\link{as.Surv}}).
##' This should follow a straight line through the origin with a slope of \eqn{1}.
##' \cr \cr
##' To check if the coefficients are close to their true values,
##' we can compute the Nelson-Aalen estimator of the
##' cumulative hazard rate of the \eqn{r_i}{r[i]}'s. A plot of this
##' estimator against \eqn{r_i}{r[i]} should be a straight line
##' through the origin with a slope of \eqn{1}.
##' \cr \cr
##' Continuous predictors are split into quantiles to facilitate the following plots:
##' \itemize{
##'   \item Plots of time vs. cumulative hazard, per predictor,
##'         should be a constant multiples of a baseline hazard, i.e. parallel.
##'   \item Plots of time vs. difference in log hazards, per predictor,
##'         should be constant over time i.e. parallel.
##'         \cr
##'         The difference should be close to \eqn{0}.
##'   \item Plots of hazard vs. reference group, per predictor,
##'         should be linear with a slope of 45 degrees.
##' }
##' 
##' \bold{Discretizing a continuous variable}
##' \cr \cr
##' These methods work by stratifying a covariate \eqn{K} into
##' \eqn{q} disjoint quantiles or groups \eqn{g}.
##' A stratified \code{coxph} model is fitted to these
##' quantiles and one is selected as a reference.
##' \cr
##' The \bold{cumulative hazard} \eqn{\hat{H}_g(t)}{H[g](t)} is
##' plotted for each group \eqn{g}.
##' These should be a constant multiple of the reference stratum
##' \eqn{\hat{H}_1(t)}{H[1](t)} over time.
##' \cr \cr
##' A simpler way to compare these is to plot the
##' \bold{differences in log cumulative hazard}, that is:
##' \deqn{\log \hat{H}_g(t) - \log \hat{H}_1(t), \quad g = 2,...,q}{
##'        log(H[1](t)) - log(H[g](t)) for g = 2 to q}
##' Each curve should be horizontal and constant over time.
##' \cr
##' Curves above zero indicate
##' an increase in hazard in group \eqn{g} vs. the reference at that time.
##' \cr \cr
##' \bold{Andersen plots} show \eqn{\log \hat{H}_1(t)}{log(H[1](t))}
##' vs. \eqn{\log \hat{H}_g(t), \quad g = 2,...,q}{
##'          log(H[g](t)) , g = 2 to K}.
##' \cr
##' If proportional hazards are present, these should be straight lines through
##' the origin.
##' \cr
##' If the curve is convex (towards the upper left of the plot) this shows that
##' \eqn{ \hat{H}_g(t) \div \hat{H}_1(t)}{
##'        H[g](t) / H[1](t)}
##' is an increasing function of \eqn{t}{time}.
##' Thus if convex, the hazard rate in group \eqn{g} is
##' increased vs. the reference, group \eqn{1}.
##' \cr \cr
##' A model with \bold{time-dependent coefficients} should not vary from one
##' without such coefficients if the assumption of proportional-hazards is met. That is, making
##' the coefficient a function of time,
##' \eqn{\hat{\beta}_k \rightarrow f(\hat{\beta}_k, t)}{
##'      B[k] --> f(B[k], t)}
##' and plotting this against time \eqn{t} should give a horizontal line.
##' \cr
##' To test this we plot the \emph{scaled Schoenfeld residuals} against time.
##' These are
##' \deqn{s^*_i = V^{-1}(\hat{\beta}, t_i)s_i}{
##'       s*[i] = INVERSE(Var(B, t)) . s[i]}
##' 
##' \subsection{A note on generating Schoenfeld residuals}{
##'   These are based on the contribution of each observation to
##'   the derivative of the log partial likelihood.
##'   \cr
##'   They are defined for each time where an event occurs and have a value
##'   for each coefficient \eqn{\hat{\beta}_k}{B[k]}. They are given by:
##'   \deqn{s_{ik} = X_{ik} - \bar{x}_k, \quad i=1,...,n \quad k=1,...,p}{
##'         s[i,k] = X[i,k] - x[k] for i=1 to n, k=1 to p}
##'   Here, \eqn{\bar{x}_k}{x[k]} is the mean of those still at risk for covariate \eqn{k}.
##'   \cr
##'   This is a weighted mean of the values of \eqn{X_k}{X[k]}
##'   (for cofficient \eqn{k}):
##'   \deqn{\bar{x}_k = \frac{\sum W X_i}{\sum W}}{
##'         x[k] = SUM W.X[i] / SUM W}
##'   and the weights are:
##'   \deqn{W = \exp \hat{\beta} X_i}{
##'         W = exp( B.X[i] )}
##'   where \eqn{X_i}{X[i]} refers to those still at risk at time \eqn{t_i}{t[i]}.
##'   \cr
##'   Now the inverse of the variance of \eqn{s_{ik}}{s[ik]} is approximately:
##'   \deqn{V^{-1} \approx Y\hat{V(\beta)}}{
##'         INVERSE(Var) = Y.Var(B)}
##'   where \eqn{Y} is the number of events and \eqn{\hat{V(\beta)}}{Var(B)} is the covariance
##'   matrix of the estimated coefficients.
##' }
##' 
##' Given the above
##' \deqn{E(s^*_i) + \hat{\beta}_k \approx \hat{\beta}_k(t_i)}{
##'       E(s*[i]) + B[k] = B(k)(t[i])}
##' so that a plot of time vs. \eqn{s^*_i} should be horizontal.
##' \cr \cr
##' \bold{TESTS OF LINEARITY FOR CONTINUOUS VARIABLES}
##' \cr \cr
##' The \bold{martingale} residual is used to help determine the best
##' functional
##' form of a covariate in a \code{coxph} model.
##' As above, the Cox model assumes that the hazard function satisfies:
##' \deqn{H_{i}(t) = H_0(t) \exp{X_i\hat{\beta}}}{
##'        H[i](t) = H[0](t) exp (X[i].B)}
##' That is, for a continuous variable,
##' a unit increase in the variable produces the same
##' change in risk across the value of the variable.
##' (E.g. an increase in age of 5 years leads to the same
##' change in hazard, no matter what the increase is from or to).
##' \cr
##' To verify this is the case, a null model is fitted
##' (i.e no coefficients, similar to intercept-only model in linear regression).
##' Martingale residuals are calcuated for this.
##' \cr
##' Plots of these residuals against the values of each of the
##' predictors in the model are shown.
##' If the correct model for covariate \eqn{k}
##' is based on a smooth function \eqn{f()}, i.e.
##' \eqn{ \exp(f(X_k)\hat{\beta_k)}}{
##'        exp(f(X[k]).B[k])}
##' then the following should hold:
##' \deqn{E(M_i | X_{ik}=X_k) \approx c.f(X_k)}{
##'       E( M[i] | X[ik]=X[k] ) = c.f(X[k])}
##' Where \eqn{M_i}{M[i]} is the martingale residual and
##' the constant \eqn{c} depends
##' on the amount of \bold{c}ensoring and is
##' roughly independent of \eqn{X_k}{X[k]}.
##' \cr
##' A \code{lowess} smoothed line is added to the plot.
##' This should be approximately linear if the
##' assumption of proportional hazards is met.
##' If the plot shows a sharp threshold, a discretised version of the covariate
##' may be preferable.
##' \cr \cr
##' \bold{Poisson} regression models \emph{are also} proportional hazards
##' models. The Cox model may thus be rewritten in Poisson form to allow
##' for application of residual methods applicable to Poisson regression.
##' \cr
##' The Cox model can be written as:
##' \deqn{H_i(t) = \exp(f(x) \hat{\beta} H_0(t))}{
##'       H[i](t) = \exp( f(x)B.H[0](t) )}
##' And the standard Poisson model can be written as:
##' \deqn{ E(Y_i|X) = \exp X_iT}{
##'        E( Y[i] | X ) = exp( X[i].T )}
##' where \eqn{Y_i}{Y[i]} are the observed events and
##' \eqn{T} is the observation time (often referred to as \eqn{\theta}{theta}).
##' This is thus:
##' \deqn{E(Y_i|X) = \exp( (f(X_i) \hat{\beta}) \int_0 Y_i(t) H_0(t) dt)}{
##'       E( Y[i] | X ) = exp( (f(X[i]).B) . INTEGRAL Y[i](t).H[0](t) dt )}
##' Where \eqn{T = \int_0 Y_i(t) H_0(t) dt}{T = INTEGRAL Y[i](t).H[0](t) dt}.
##' Once expressed as a Poisson model, this can be analysed using tools
##' available for generalized additive models (\code{gam}).
##' \cr
##' To do this the \code{coxph} model is refit with \code{gam}. The
##' outcome (left-hand side) is those times in which an event occurred.
##' The predictors are the same. For continuous terms, an attempt is
##' made to fit a non-linear function to improve the fit, with a default
##' of \eqn{4} degrees of freedom (\code{degfP=4}).
##' \cr
##' The Poisson model fit with \code{gam} has an additional \code{offset}
##' term. This is
##' \deqn{\texttt{offset} = \log[ \exp(-X_i\hat{\beta}) \exp(X_i\hat{\beta}H_0(t_i))]
##'              = H_0(t_i)
##'              \approx  \int_0 Y_i(t)H_0(t) dt}{
##' offset = log( exp(-X[i].B) . exp(X[i].B.H[0](t)) ) = INTEGRAL Y[i](t)H[0](t) dt }
##' See \code{?predict.coxph} for details.
##' Plots and anova are generated (see \code{?anova.gam} for details).
##' Plots show the residuals by the values of the predictor.
##' Ideally these should be horizontal with the intercept at zero.
##' \cr \cr
##' \bold{Regression splines} may be used to replace continuous terms
##' directly in the \code{coxph} function.
##' These are fit by connecting number of knots with locally fitting curves.
##' The degrees of freedom (by default \code{degfRS=4}) is the number of knots plus one.
##' \code{degfRS} \eqn{-1} dummy variables are generated to try to improve the
##' fit of the variable. Plots of the original variable vs. the fitted splines
##' should be linear. The function uses B-splines.
##' \cr \cr
##' \bold{Penalized smoothing splines} are an alternative to
##' regression splines which, for small degrees of freedom, have better
##' properties regarding their locality of influence.
##' They are chosen to minimize \eqn{\beta}{B} for the basis functions, in:
##' \deqn{ \theta \sum_{i=1}^n [y_i -f(x_i, \beta)]^2 +
##'          (1-\theta) \int [f''(x,\beta)]^2 dx }{
##'        T . SUM [y - f(x, B)]^2 + (1-T) . INTEGRAL [f''(x, B)]^2 dx} 
##' Here the first term is the residual sum of squares and the
##' second is the integral of the second derivative of the function \eqn{f}
##' with respect to \eqn{x}.
##' \cr
##' For a straight line \eqn{f''(x)=0} and the term will increase in
##' proportion to the degree of curvature.
##' \eqn{\theta}{T} is a \bold{t}uning parameter based on the degrees of freedom
##' (by default \code{degfSS=4}).
##' As \eqn{\theta \rightarrow 0}{T approaches zero}
##' (\eqn{2} degrees of freedom, including intercept),
##' the solution converges to the least-squares line.
##' As \eqn{\theta \rightarrow 1}{T approaches one}, (\eqn{n} degrees of freedom),
##' the solution approaches a curve that passes through each point.
##' Plots of the fitted splines vs. the original variable
##' should be linear.
##' \cr \cr
##' \bold{TESTS OF INFLUENCE}
##' \cr \cr
##' The simplest measure of influence is the \bold{jackknife} value
##' \deqn{J_i=\hat{\beta} - \hat{\beta}_{-i}}{
##'       J[i] = B - B[-i]}
##' where \eqn{\beta_{-i}}{B[-i]} is the result of a fit that includes all
##' observations except \eqn{i}. This can be computed as
##' \itemize{
##'  \item Converge to \eqn{\hat{\beta}}{B} as usual e.g. via Newton-Raphson.
##'  \item Delete observation \eqn{i}.
##'  \item Perform one additional iteration.
##' }
##' This may be expressed as:
##' \deqn{\delta \beta = 1'(U \chi^{-1}) = 1'D}{
##'       change in B = 1'(U.INVERSE(Var(B))) = 1'D} 
##' Here \eqn{D}, the matrix of \bold{dfbeta residuals}, is comprised of
##' the score residuals \eqn{U} scaled by the variance of \eqn{\beta}{B},
##' \eqn{var(B)}{Var(B)}. Each row of \eqn{D} is the change in \eqn{\hat{\beta}}{B} if
##' observation \eqn{i} is removed.
##' \cr \cr
##' Caution - for plots to verify proportional hazards:
##' the variance of the curves is not constant over time.
##' Increasing departures from model assumptions are likely to be found
##' as time increases. 
##' \code{package:surv2sample} may need to be installed from source
##' to allow one of the examples to run.
##'
##' @references Examples are from
##' \bold{K&M} Example 11.1 - 11.7, pg 355--66.
##' @references Last example is from:
##' Therneau T, Grambsch P 2000.
##' \emph{Modeling Survival Data}, 1st edition.
##' New York: Springer.
##' Section 5.1.2, pg 91.
##' \href{http://dx.doi.org/10.1007/978-1-4757-3294-8}{Springer (paywall)}
##' @references
##' Andersen PK, Borgan O, Gill R, Keiding N 1982.
##' Linear Nonparametric Tests for Comparison of Counting Processes, with
##' Applications to Censored Survival Data, Correspondent Paper.
##' \emph{International Statistical Review} \bold{50}(3):219--44.
##' \href{http://www.jstor.org/stable/1402489}{JSTOR}
##'
##' @examples
##' \dontrun{
##' ### running these examples with toPdf=FALSE will 
##' ### open too many devices to be compatible with R CMD check
##' ### results from these examples can be found in the package source
##' ### under /survMisc/inst/doc/
##' ### 
##' ### for log-log plot
##' if(require(devtools)){
##' ### this package is now archived, so need to install from url
##'  install_url("http://cran.r-project.org/src/contrib/Archive/surv2sample/surv2sample_0.1-2.tar.gz")
##'      library(surv2sample)
##'      data(gastric, package="surv2sample")
##'      dx(coxph(Surv(time/365, event) ~ treatment, data=gastric), file="gasDx.pdf")
##'  }
##' data(bmt, package="KMsurv")
##' bmt <- within(bmt, {
##' z1 <- z1 -28
##' z2 <- z2- 28
##' z3 <- z1*z2
##' z4 <- as.double(group == 2)
##' z5 <- as.double(group == 3)
##' z6 <- z8
##' z7 <- (z7 / 30) - 9
##' z8 <- z10
##' })
##' c1 <- coxph(Surv(t2, d3) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8,
##'          method="breslow", data=bmt)
##' dx(c1, file="bmtDx.pdf")
##' ###
##' data(alloauto, package="KMsurv")
##' c1 <- coxph(Surv(time,delta) ~ factor(type),
##'             method="breslow", data=alloauto)
##' dx(c1, file="alloDx.pdf")
##' ### GAM model. Therneau 5.3
##' data(pbc, package="survival")
##' w1 <- which(is.na(pbc$protime))
##' pbc <- pbc[-w1, ]
##' c1 <- coxph(Surv(time, status==2) ~ age + edema + bili + protime + albumin,
##'             data=pbc, method="breslow")
##' dx(c1, file="pbcDx.pdf")
##' ### Time dependent covariate. Therneau 6.3
##' data(veteran, package="survival")
##' veteran$celltype <- relevel(veteran$celltype, ref="adeno")
##' c1 <- coxph(Surv(time, status) ~ trt * celltype + karno + diagtime + log(age) + prior,
##'             data=veteran[-1, ])
##' dx(c1, what="ph", file="vetDx.pdf")
##' }
##' ### simple example which doesn't take up too many devices
##' c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
##' dx(c1)
##' 
dx.coxph <- function(x,
                     ...,
                     what=c("all", "ph", "lin", "inf"),
                     toPdf=TRUE,
                     file="dxPlots.pdf",
                     maxStrata=5,
                     defCont=2,
                     noQuantiles=3,
                     maxFact=4,
                     identify=FALSE,
                     degfP=3,
                     degfRS=4,
                     degfSS=4,
                     timeTrans=c("km", "log", "rank", "identity"),
                     ties
                     ){
    stopifnot(inherits(x, "coxph"))
    if (toPdf && identify) stop("Cannot use 'identify' (interactive) when writing to file")
    makeTitle <- function(main1) graphics::mtext(main1, line=0.3, outer=TRUE)
### for R CMD check
    e <- E <- n <- offT1 <- d1 <- NULL
### 
    if (toPdf) grDevices::pdf(file)        
###
    if(missing(ties)) ties <- x$method
    stopifnot(ties %in% c("breslow", "efron", "exact"))
    timeTrans <- match.arg(timeTrans)
    what <- match.arg(what)
    res1 <- vector(mode="list")
### some duplication below, not v. efficient...
### 
### no. events (take from 2nd column of Surv object)
    e1 <- model.frame(x)[, 1][, "status"]
    t1 <- model.frame(x)[, 1][, "time"]
    dt2 <- tne(x)[, c("t", "n", "e"), with=FALSE]
###
    if(what=="all" | what=="ph") {
###
###----------------------------------------
###----------------------------------------
### Tests of proportional hazards
###----------------------------------------
###----------------------------------------
### 
        main2 <- paste0("Complete model:", "\n", paste0(deparse(x$formula), collapse=" "))
### ### remove whitespace
        main2 <- gsub("  *", " ", main2)
### 
###----------------------------------------
### log(log) plot 
###----------------------------------------
###
### ### ensure not too many strata for plot
        if (length(unique(tne(x)[, s])) < maxStrata){
            main1 <- "Time vs. -log-log survival.
Kaplan-Meier curves should be parallel if
proportional-hazards assumption is correct"
### ### ### change coxph to survfit
            c2 <- paste(deparse(x$call), collapse="")
### ### ### get name of data from call
            d1 <- x$call$data
            storage.mode(d1) <- "character"
            f1 <- deparse(x$formula)
            f2 <- paste("survfit(", f1, ", data=", d1, ")")
            s1 <- eval(parse(text=f2))
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            plot(s1$time,
                 -log(-log(s1$surv)),
                 xlab="Time",
                 ylab="-log(-log(Survival))",
                 main=main2,
                 ...)
            makeTitle(main1=main1)
        }
###
###----------------------------------------
### Expected = Cox-Snell residuals
###----------------------------------------
###
        dt2[, E:= predict(x, type="expected")]
### or from survival:::residuals.coxph()
### coxSnellRes1 <- dt2[, e] - residuals(x, type="martingale")
### stopifnot(predict(x, type="expected")==coxSnellRes1)
###
###----------------------------------------
### compare to unit exponential distribution
###----------------------------------------
###         
        main1 <- "Quantile-quantile plot.
Unit exponential distribution vs. expected events (Cox-Snell residuals).
Should follow line through origin at 45 degrees (blue) if well fit."
        graphics::par(oma=c(0, 0, 4, 0))
        set.seed(1)
        stats::qqplot(y=dt2[, E],
                      x=rexp(nrow(dt2)),
                      ylab="Expected events (Cox-Snell residuals)",
                      xlab="Unit exponential distribution",
                      main=main2,
                      ...)
### 
### alternatively can show axes on log scale for better display
### add the following to qqplot above
### log="xy", axes=FALSE
### then 
### axis(1, at=round(q1$x, 2)[seq(from=1, to=length(q1$x), length.out=11)])
### axis(2, at=round(q1$y, 2)[seq(from=1,to=length(q1$y), length.out=11)])
### 
        stats::qqline(y=dt2[, E],
                      distribution=qexp,
                      col="blue", lty=2)
        makeTitle(main1=main1)
###
###----------------------------------------
### observed vs. expected hazard
###----------------------------------------
### 
        main1 <- "Observed vs. expected hazard.
Should follow line through origin at 45 degrees (blue) if well fit."
        h1 <- sf(e=dt2[, sum(e), by=t]$V1,
                 n=dt2[, max(n), by=t]$V1, what="h")
        h2 <- sf(e=tne(as.Surv(dt2[, t], dt2[, E]))$e,
                 n=tne(as.Surv(dt2[, t], dt2[, E]))$n,
                 what="h")
        if(!toPdf) dev.new()
        graphics::par(oma=c(0, 0, 4, 0))
        graphics::plot(h1, h2,
                       xlab="Observed hazard",
                       ylab="Expected hazard",
                       main=main2,
                       ...)
        abline(0, 1, col="blue", lty=2)
        makeTitle(main1=main1)
### 
###----------------------------------------
### cumulative hazard based on expected values
###----------------------------------------
### 
        main1 <- "Expected events vs. hazard based on sorted expected events
or Cox-Snell residuals vs. cumulative hazard of these residuals.
Should follow line through origin at 45 degrees (blue) if well fit."
### 
        s1 <- survival::survfit(Surv(dt2[, E], dt2[, e]) ~ 1)
        h3 <- sf(e=s1$n.event, n=s1$n.risk, what="h")
### ### similarly
### ### stopifnot(dt2[, sum(e), by=E][order(E)]$V1==s1$n.event)
        if(!toPdf) dev.new()
        graphics::par(oma=c(0, 0, 4, 0))
        plot(x=s1$time,
             y=h3,
             type='s',
             xlab="Expected events (sorted), or Cox-Snell residuals",
             ylab="Hazard = cumulative sum of [events/ no. at risk] (sorted by expected events)",
             main=main2,
             ...)
        graphics::abline(0, 1, col="blue", lty=2)
        makeTitle(main1=main1)
### 
###----------------------------------------
### Cox-Snell stratified by each co-variate
###----------------------------------------
###
        main1 <- "Time vs. hazard, per predictor.
If hazards proportional then curves should be constant multiples of a baseline.
Reference (black) line is 45 degrees."
        for (i in 2:length(model.frame(x))){
            n1 <- names(model.frame(x))[i]
            x1 <- model.frame(x)[, i]
### ### ### split if continous (into quantiles)
            if (length(unique(x1)) > defCont && !(is.factor(x1))){
                x1 <- Hmisc::cut2(x1, g=noQuantiles)
            }
            x1 <- as.factor(x1)
            dt2[, x1 := as.factor(x1)]
            levX1 <- levels(x1)
### ### ### refit as stratified model
            f1 <- stats::as.formula("Surv(t1, e1) ~ strata(x1)")
            c2 <- survival::coxph(f1, ties=ties)
            s1 <- survival::survfit(c2, type="aalen")
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            plot(s1, fun='cumhaz',
                 ylab="Cumulative hazard",
                 xlab="Time",
                 main=paste0("Predictor: ", n1),
                 col=rainbow(length(levX1), start=0.5, end=1),
                 ...)
### ### ### add legend
            graphics::legend("bottomright", levX1,
                             title = paste0("Levels of: ", n1),
                             lty=1,
                             col=rainbow(length(levX1), start=0.5, end=1),
                             ...)
### ### ### which value of cumulative hazard closest to 0.25
### ### ### used below to plot line            
            w1 <- which(abs(s1$cumhaz - 0.25) == min(abs(s1$cumhaz - 0.25)))
### ### ### 45 degree line
            graphics::abline(0, 1 / (4 * s1$time[w1]), lty=2)
            makeTitle(main1=main1)
        }
###
###----------------------------------------
### differences in log-hazard curves
### note repeats code above which is inefficient
### but this is easier to read...
###----------------------------------------
###
        main1 <- "Time vs. difference in log hazards, per predictor.
Should be constant over time.
If >0 (black line) shows survival advantage for reference group."
### store merged frames and names
        m1l <- ns1l <- vector("list", (ncol(model.frame(x))-1))
        for (i in 2:length(model.frame(x))){
            n1 <- names(model.frame(x))[i]
            mf1 <- model.frame(x)[, i]
            if (is.factor(mf1) && length(unique(mf1)) > maxFact) next
            if (length(unique(mf1)) > defCont & !(is.factor(mf1))){
                mf1 <- Hmisc::cut2(as.numeric(mf1), g=noQuantiles)
            }
### ### ### put original time + events into data.frame
            df1 <- data.frame(t1=unclass(model.frame(x)[[1]])[, 1],
                              e1=unclass(model.frame(x)[[1]])[, 2],
                              mf1=mf1)
            df1$mf1 <- as.factor(df1$mf1)
            ns1 <- levels(df1$mf1)
### ### ### refit as stratified model
            f1 <- stats::as.formula("Surv(t1, e1) ~ strata(mf1)")
            c2 <- survival::coxph(f1, data=df1, ties=ties)
            s1 <- survfit(c2, type="aalen")
### ### ### get baseline hazard (= -log survival)
            s2 <- data.frame(t=s1$time, -log(s1$surv))
### ### ### names of strata
            ns1 <- paste0(n1, "=", ns1 )
            s2$strat <- unlist(sapply(1:length(s1$strata), function(i) rep(ns1[i], s1$strata[i]) ) )
### ### ### split by strata
            s3 <- by(s2[, -3], s2[, 3], FUN=identity)
### ### ### rename 2nd column to avoid duplicate column names
### ### ### needed if merging >3 data.frames
            for(j in seq_along(s3)) {
                colnames(s3[[j]])[2] <- ns1[j]
            }
### ### ### merge list based on time
            m1 <- Reduce(function(...)
                         merge(..., by="t",
                               all=TRUE), s3)
### ### ### replace NAs in first row with zeros
            m1[1, is.na(m1[1, ])] <- 0
### ### ### last observation carried forward
            m1 <- zoo::na.locf(m1)
### ### ### get differences between reference and other strata
            diffs1 <- sapply(3:ncol(m1), function (i) -log(m1[, 2]) + log(m1[, i]))
### ### ### one line per difference
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            graphics::matplot(m1$t, diffs1,
                              type="s",
                              lty=1,
                              xlab="Time",
                              ylab="Difference in log cumulative hazard",
                              col=rainbow(length(ns1), start=0.7, end=1),
                              main=paste0("Reference: ", gsub("=", " = " , ns1[1])),
                              ...)
### ### ### add legend
            l1 <- sapply(1:ncol(diffs1),
                         function(i) paste0(ns1[i+1], "  -  ", ns1[1])
                         )
            tit1 <- "Group - reference group: "
            graphics::legend("bottomright", legend=l1,
                             lty=1,
                             title=tit1,
                             cex=0.75,
                             col=rainbow(length(ns1), start=0.7, end=1)
                             )
            graphics::abline(0, 0, lty=2)
            makeTitle(main1=main1)
### ### ### store results in list for Anderson plot
            m1l[[i-1]] <- m1
            ns1l[[i-1]] <- ns1
        }
###
###----------------------------------------
### Andersen plots
###----------------------------------------
###
        main1 <- "Cumulative hazard vs. reference group. Should be linear plot through origin.
If convex (towards top left) shows ratio of hazards is increasing over time.
Reference line (black) is at 45 degrees."
        for (i in 2:length(model.frame(x))){
            n1 <- names(model.frame(x))[i]
            m1 <- m1l[[i-1]]
            ns1 <- ns1l[[i-1]]
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            matplot(m1[, 2], m1[, -c(1,2)],
                    lty=1,
                    xlim=c(0,1),
                    xlab=paste0("Cumulative hazard (reference): \n", ns1[1]),
                    ylim=c(0,1),
                    ylab="Cumulative hazard (see legend)",
                    type="s",
                    col=rainbow(length(ns1)-1, start=0.7, end=1),
                    main=paste0("Predictor: ", n1),
                    ...)
            abline(0, 1, lty=2)
            graphics::legend("topleft", legend=ns1[-1],
                             lty=1,
                             col=rainbow(length(ns1[-1]), start=0.7, end=1)
                             )
            makeTitle(main1=main1)
        }
    }
###
    if(what=="all" | what=="lin"){        
###----------------------------------------
###----------------------------------------
### Tests of linearity of predictors
###----------------------------------------
###----------------------------------------
###
###----------------------------------------
### Poisson (GAM) plots
###----------------------------------------
###
### ### adapted from code for lm()
        mf <- x$call
        m <- match(c("formula", "data"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
### ### model terms
        mt <- attr(mf, "terms")
        stopifnot(is.empty.model(mt)==FALSE)
### ### if at least some terms, exclude intercept
        if (length( attr(mt, "term.labels") ) > 0) {
            attr(mt, "intercept") <- 0
        }
### ### get model matrix
        dt1 <- data.table(model.matrix(mt, mf, contrasts))
### ### need unlist to prevent binding later of names later
        n1 <- unlist(list(names(dt1)))
### ### check which are factors
        fact1 <- names(attr(terms(x), "dataClasses")=="factor")[
                               attr(terms(x), "dataClasses")=="factor"]
        if(length(fact1)){
### ### replace factor levels in model.matrix with original factor
### grep1 <- paste0("^", fact1, ".*")
            w1 <- which(grepl(fact1, n1, fixed=TRUE))
            dt1[, w1:=NULL, with=FALSE] 
            if (ncol(dt1)) dt1[, eval(quote(fact1)) := eval(x$call$data)[, fact1], with=FALSE]
        }
### ### check which have at least 3 values/ levels and are not factors 
        gr3 <- sapply(dt1, function(i) length(unique(i)) > 3 && !is.factor(i))
        if (any(gr3)) {
            main1 <- "Poisson approach to determine non-linearity.
Predictor (with smoothing spline) vs. residuals from GAM plot.
If linear should be horizontal line with intercept=0 (blue)."
            n1 <- unlist(list(names(dt1)))
### ### ### replace names with x1, x2...xn to build formula
### ### ### for variable names like e.g. 'log(x1)' and 'x1:x2'
            n2 <- paste0("x", seq(ncol(dt1)))
            data.table::setnames(dt1, n2)
            p1 <- paste("s(", n2[gr3], ", ", degfP, ")", sep="", collapse=" + ")
            p2 <- paste(n2[!gr3], sep="", collapse=" + ")
            p1 <- paste(p1, p2, sep=" + ")
### ### ### add offset
            dt1[, offT1 :=  predict(x, type="expected") * exp(-predict(x, type="lp"))]
### ### ### add column for events
            dt1[, e1 := e1]
### ### ### make formula
            p3 <- paste("e1 ~ ", p1, "+ offset(log(offT1))")
### ### ### GAM model
            f2 <- paste("gam(", p3, ", data=dt1, family=poisson)")
            g1 <- eval(parse(text=f2))
### ### ### from gam:::anova.gam
            a1 <- anova(g1)
### ### ### get continuous/ smoothed terms
            te1 <- rownames(a1)[grep("s()", rownames(a1))]
            yn1 <- mapply(gsub, pattern=n2[gr3], replacement=n1[gr3], x=te1)
### ### ### from gam:::preplot.gam
            pg1 <- preplot(g1, rug=TRUE, se=TRUE, terms=te1)
            for (i in seq_along(pg1)){
                if(!toPdf) dev.new()
                graphics::par(oma=c(0, 0, 4, 0))
                main2 <- paste0("Predictor: ", yn1[i]) 
### ### ### ### scatterplot is invisible
                with(pg1[[i]], graphics::plot(x, y, col="white",
                                              xlim=c(min(x), max(x)),
                                              ylim=c(min(y - se.y), max(y + se.y)),
                                              xlab=yn1[i],
                                              ylab=paste0("Partial residual for ", yn1[i]),
                                              main=main2,
                                              ...))
### ### ### ### loess smoother
### ### ### ### span=0.1 is arbitrary
### ### ### ### but closely approximates apprearance of plot.gam()
                lo1 <- with(pg1[[i]], loess(y ~ x))
### ### ### ### additional interpolated points
                xl1 <- with(pg1[[i]], seq(min(x), max(x), (max(x) - min(x))/100))
                graphics::lines(xl1, predict(lo1, xl1), col='black', lwd=2)
### ### ### ### add std errors
                lo1 <- with(pg1[[i]], loess(y + se.y ~ x))
                graphics::lines(xl1, predict(lo1, xl1), col='black', lty=2)
                lo1 <- with(pg1[[i]], loess(y - se.y ~ x))
                graphics::lines(xl1, predict(lo1, xl1), col='black', lty=2)
### ### ### ### add rug to x axis
                with(pg1[[i]], rug(x))
                graphics::abline(0, 0, col="blue")
                makeTitle(main1=main1)
            }
            res1[["Poisson"]] <- a1
        } else {
            res1[["Poisson"]] <- "No predictors with >3 unique values"
        }
### 
###----------------------------------------
### regression splines
###----------------------------------------
###
        if(ncol(dt1)) suppressWarnings(dt1[, offT1 := NULL])
        if(ncol(dt1)) suppressWarnings(dt1[, e1 := NULL])
        grDegfRS1 <- apply(dt1, 2, function(i) length(unique(i)) > degfRS)
        if (any(grDegfRS1)){
            main1 <- "Regression splines approach to determine non-linearity.
Variable vs. partial residual for smoothed variable.
If linear should be horizontal line with intercept=0 (blue)"
### ### ### proper names (for plot); i.e. not x1, x2, ...
            xna1 <- n1[grDegfRS1]
### ### ### names used in formula
            xna2 <- paste("ns(", n1[grDegfRS1],
                          ", df=", degfRS, ")", sep="")
### ### ### for formula
            p1 <- paste("ns(", n2[grDegfRS1],
                        ", df=", degfRS, ")", sep="", collapse=" + ")
            p2 <- paste(n2[!grDegfRS1], sep="", collapse=" + ")
            p3 <- ifelse(p2=="", p1, paste(p1, p2, sep=" + "))
            p3 <- paste("Surv(t1, e1) ~", p3)
            dt1[, t1 := t1]
            dt1[, e1 := e1]
            f2 <- paste("coxph(", p3, ", data=dt1)")
            c2 <- eval(parse(text=f2), envir=environment())
### ### ### from survival:::predict.coxph
            p1 <- predict(c2, type="terms")
            for (i in 1:length(which(grDegfRS1==TRUE))){
                if(!toPdf) dev.new()
                graphics::par(oma=c(0, 0, 4, 0))
                main2 <- paste0("Predictor: ", xna1[i], " with ", degfRS, " df")
                plotTerm(c2, term=i,
                         xlab=xna2[i],
                         ylab=paste0("Partial residual for: ", xna1[i]),
                         data=dt1,
                         main=main2)
                graphics::abline(0, 0, col="blue")
                makeTitle(main1=main1)
            }
        }
###
###----------------------------------------
### smoothing splines
###----------------------------------------
###
        if(ncol(dt1)) suppressWarnings(dt1[, t1 := NULL])
        if(ncol(dt1)) suppressWarnings(dt1[, e1 := NULL])
        grDegfSS1 <- apply(dt1, 2, function(i) length(unique(i)) > degfSS)
        if (any(grDegfSS1)){     
            main1 <- "Smoothing splines approach to determine non-linearity.
Variable vs. partial residual for smoothed variable.
If linear should be horizontal line with intercept=0 (blue)"
### ### ### proper names (for plot)
            xna1 <- n1[grDegfSS1]
            xna2 <- paste("pspline(", n1[grDegfSS1],
                          ", df=", degfRS, ")", sep="")
            p1 <- paste("pspline(", n2[grDegfSS1],
                        ", df=", degfSS, ")", sep="", collapse=" + ")
            p2 <- paste(n2[!grDegfSS1], sep="", collapse=" + ")
            p3 <- ifelse(p2=="", p1, paste(p1, p2, sep=" + "))
            p3 <- paste("Surv(t1, e1) ~", p3)
            dt1[, t1 := t1]
            dt1[, e1 := e1]
            f2 <- paste("coxph(", p3, ", data=dt1)")
            c2 <- eval(parse(text=f2))
            for (i in 1:length(which(grDegfSS1==TRUE))){
                if(!toPdf) dev.new()
                graphics::par(oma=c(0, 0, 4, 0))
                main2 <- paste0("Predictor: ", xna1[i], " with ", degfSS, " df")
                plotTerm(c2, i,
                         xlab=xna2[i],
                         ylab=paste0("Partial residual for ", xna1[i]),
                         data=dt1,
                         main=main2)
                graphics::abline(0, 0, col="blue")
                makeTitle(main1=main1)
            }
        }
###
###----------------------------------------
### time-dependent coviate
###----------------------------------------
###
        z1 <- survival::cox.zph(x, transform=timeTrans)
        main1 <- paste0("Time vs. scaled Schoenfeld residuals, with smoothed spline (black).
If <0 (red line), indicates a protective effect.
Regression line (blue) should be horizontal if model well fit.")
        max1 <- ifelse("GLOBAL" %in% rownames(z1$table),
                       nrow(z1$table)-1,
                       nrow(z1$table))
        for (i in seq(max1)){
            main2 <- paste0("Predictor: ", rownames(z1$table)[i], "
Time transform: ", timeTrans)
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            graphics::plot(z1[i], main=main2)
            l1 <- stats::lm(z1$y[, 1] ~ z1$x)
            graphics::abline(a=l1$coef[1], b=l1$coef[2], col="blue")
            graphics::abline(a=0, b=0, col="red")
            makeTitle(main1=main1)
        }
        res1[[paste0("tt, ", timeTrans)]] <- z1
    }
###
    if(what=="all" | what=="inf"){
###
###----------------------------------------
###----------------------------------------
### Tests of influence
###----------------------------------------
###----------------------------------------
###
###----------------------------------------
### DFbeta
###----------------------------------------
###
        main1 <- "Coefficient vs. jacknife influence.
Change in coefficient if this observation dropped.
Outliers may need to be re-examined"
        dfb1 <- as.matrix(residuals(x, type="dfbeta"))
        colnames(dfb1) <- names(x$coefficients)
        mm1 <- model.matrix(x)
        for (i in 1:length(x$coefficients)){
            n1 <- names(x$coefficients)[i]
            main2 <- paste0("Coefficient: ", n1)
            if (!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            x1 <- mm1[, n1]
            if (length(unique(x1)) > defCont){
                graphics::plot(x1, dfb1[, i],
                               xlab=n1,
                               ylab="Influence (dfBeta)",
                               main=main2)
            } else {
                graphics::boxplot(dfb1[, i] ~ x1,
                                  col="lightgrey",
                                  xlab=n1,
                                  ylab="Influence (dfBeta)",
                                  main=main2)
            }
            graphics::abline(a=0, b=0)
            makeTitle(main1=main1)
            if(identify) identify(x1, dfb1[, i], labels=rownames(dt1))
        }
###
###
### 
        main1 <- "Coefficient vs. jacknife influence scaled by standard error of coefficients.
Change in coefficient if this observation dropped.
Outliers may need to be re-examined"
### "dfbetas" gives scaled (vs. "dfbeta" above)
        dfb1 <- as.matrix(residuals(x, type="dfbetas"))
        colnames(dfb1) <- names(x$coefficients)
        for (i in 1:length(x$coefficients)){
            n1 <- names(x$coefficients)[i]
            main2 <- paste0("Coefficient: ", n1)
            if (!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            x1 <- mm1[, n1]
            if (length(unique(x1)) > defCont){
                graphics::plot(x1, dfb1[, i],
                               xlab=n1,
                               ylab="Scaled influence (dfBeta)",
                               main=main2,
                               ...)
            } else {
                graphics::boxplot(dfb1[, i] ~ x1,
                                  col="lightgrey",
                                  xlab=n1,
                                  ylab="Scaled influence (dfBeta)",
                                  main=main2)
            }
            graphics::abline(a=0, b=0)
            makeTitle(main1=main1)
            if(identify) identify(x1, dfb1[, i], labels=rownames(dt1))
        }
###
###----------------------------------------
### scaled score residuals
###----------------------------------------
###
        ss1 <- residuals(x, type="score") %*% x$var
        colnames(ss1) <- names(x$coefficients)
### need to remake XS
###dt1 <- data.table(model.matrix(x))
        main1 <- "Coefficient vs. scaled score residuals.
Assesses leverage: influence of observation on a single coefficient.
Outliers may need to be re-examined"
        for (i in 1:length(x$coefficients)){
            n1 <- names(x$coefficients)[i]
            main2 <- paste0("Coefficient: ", n1)
            x1 <- mm1[, n1]
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            if (length(unique(x1)) > defCont){
                graphics::plot(x1, ss1[, i],
                               xlab=n1,
                               ylab="Leverage (scaled score residuals)",
                               type="p",
                               main=main2,
                               ...)
            } else {
                graphics::boxplot(ss1[, i] ~ x1, col="lightgrey",
                                  xlab=n1,
                                  ylab="Leverage (scaled score residuals)",
                                  main=main2)
            }
            graphics::abline(a=0, b=0)
            makeTitle(main1=main1)
            if(identify) identify(x1, ss1[, i], labels=rownames(dt1))
        }
###
###----------------------------------------
### likelihood displacement vs martingale
###----------------------------------------
###
        ld1 <- ss1 * residuals(x, type="score")
        m1 <- residuals(x, type="martingale")
        main1 <- "Martingale residuals vs. likelihood displacement residuals.
Assesses influence of observation on coefficient.
Outliers may need to be re-examined."
        for (i in 1:length(x$coefficients)){
            n1 <- names(x$coefficients)[i]
            main2 <- paste0("Coefficient: ", n1)
            if(!toPdf) dev.new()
            graphics::par(oma=c(0, 0, 4, 0))
            graphics::plot(m1, ld1[, i],
                           xlab="Martingale residuals",
                           ylab=paste0("Likelihood displacement for: ", n1),
                           type="p",
                           main=main2)
            graphics::abline(a=0, b=0)
            makeTitle(main1=main1)
            if(identify) identify(m1, ld1[, i], labels=rownames(dt1))
        }
###
###----------------------------------------
###
    }
    if(toPdf) dev.off()
    if(length(res1)) return(res1)
}
