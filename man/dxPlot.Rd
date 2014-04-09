\name{dxPlot}
\alias{dxPlot}
\alias{dxPlot.coxph}
\title{Diagnostic plots for coxph model}
\usage{
dxPlot(x, ...)

\method{dxPlot}{coxph}(x, ..., ties = "breslow", defCont = 2,
  noQuantiles = 5, noPerPage = 2, height = NULL, width = NULL)
}
\arguments{
  \item{x}{An object of class \code{coxph}}

  \item{...}{Additional arguments}

  \item{ties}{Method of handling ties when refitting model.
  Must be one of \code{breslow} or \code{efron}}

  \item{defCont}{Definition of continuous variable.
  Variables with more than \emph{n} unique values will be
  split into quantiles to facilitate graphs. (This does not
  apply to factor variables)}

  \item{noQuantiles}{No. of quantiles into which to split
  continuous variables}

  \item{noPerPage}{Number of plots per page}

  \item{width}{Width of screen (display device) in pixels.
  Set to \code{NULL} for default plot size}

  \item{height}{Height of screen (display device) in
  pixels. Set to \code{NULL} for default plot size}
}
\value{
Graphs as described above. Plotted with base graphics.
}
\description{
Diagnostic plots for coxph model
}
\details{
The \bold{Cox-Snell} residuals are used to assess the fit
of a proportional hazards model. \cr The residuals are
generated from the (non time-dependent) covariates \eqn{Z},
a matrix with one row per observation (total \eqn{n}) and
additional indicators of time \eqn{T} and status
\eqn{\delta}{D}. estimated coefficients \eqn{b}, where
\eqn{b} is a vector of length \eqn{p} (no. of predictors) :
\deqn{ r_j = \hat{H}_0(T_j) \exp ( \sum_{k=1}^p Z_{jk}b_k
), \quad j=1,...,n }{ r[j] = H0(T[j]) exp ( SUM Z[jk] b[k]
), j=1,...,n } Here \eqn{\hat{H}_0}{H0} is Breslows
estimator of the `baseline' hazard (i.e. all coefficients
are zero). \cr \cr If the coefficients are close to their
true values, then \eqn{r_j}{r[j]} should follow a
unit-exponential distribution, i.e.
\eqn{H_0(t)=t}{H0(t)=t}. \cr \cr To verify this, we
calculate the Nelson-Aalen estimator of the cumulative
hazard rate of the \eqn{r_j}{r[j]}s. A plot of this
estimator against \eqn{r_j}{r[j]} should be a straight line
through the origin with a slope of \eqn{1}. \cr \cr The
\bold{martingale} residual is used to help determine the
best functional form of a covariate in a \code{coxph}
model. The Cox model assumes that the hazard function
satisfies: \deqn{\lambda_{i}(t) = \lambda_0(t)
\exp(X_i\beta)}{ h[i](t) = h[0](t) exp (X[i]Beta)} That is,
for a continuous variable, an unit increase in the variable
produces the same change in risk across the value of the
variable. (E.g. an increase of age of 5 years leads to the
same change in hazard no matter what the increase is from
or to). \cr \cr To verify this is the case, a null model is
fitted (i.e no coefficients, similar to intercept-only
model in linear regression). Martingale residuals are
calcuated for this. Plots of these residuals against the
values of each of the predictors in the model are shown.
\cr \cr If the correct model for covariate \eqn{j} is based
on a smooth function \eqn{f()}, i.e. \eqn{
\exp(f(X_j)\beta_j)}{exp(f(X[j])Beta[j])} then the
following should hold: \deqn{E(M_i|X_{ij}=x_j) \approx
cf(x_j)}{ E( M[i]| X[i,j]=x[j] ) = cf(x[j])} Where \eqn{M}
is the martingale residual and \eqn{c} depends on the
amount of censoring and is roughly independent of
\eqn{x_j}{x[j]}. \cr \cr A \code{lowess} smoothed line is
added to the plot. This should be approximately linear if
the assumption of proportional hazards is met. If the plot
shows a threshold, a discretised version of the covariate
may be preferable. \cr \cr The assumption of
\bold{proportional hazards} can be checked in a number of
ways. These methods work by stratifying a covariate \eqn{G}
into \eqn{k} disjoint strata. A stratified \code{coxph}
model is fitted to these strata and one is selected as a
reference. \cr The \bold{cumulative hazard}
\eqn{\hat{H}_g(t)}{H[g](t)} is plotted for each stratum
\eqn{g}. These should be a constant multiple of the
reference stratum \eqn{\hat{H}_1(t)}{H[1](t)} over time.
\cr \cr Another way to compare these is to plot the
\bold{differences in log cumulative hazard}, that is:
\deqn{\log \hat{H}_g(t) - \log \hat{H}_1(t), \quad g = 2 ,
..., k}{ log(H[1](t)) - log(H[g](t)) for g = 2 to k} Each
curve should be horizontal and constant over time. Curves
above zero indicate increased hazard in the stratum \eqn{g}
vs the reference at that time. \cr \cr Finally
\bold{Andersen plots} show \deqn{\log \hat{H}_1(t) vs \log
\hat{H}_g(t), \quad g = 2,...,k}{ log(H[g](t)) vs
log(H[1](t)) , g = 2 to k} If proportional hazards are
present, these should be straight lines through the origin.
If the curve is convex this shows that \eqn{ \hat{H}_g(t)
\div \hat{H}_1(t)}{ H[g](t) / H[1](t)} is an increasing
function of \eqn{t}{time}. Thus if convex, the hazard rate
in \eqn{g} is increasing vs the reference.
}
\note{
Caution - for plots to verify proportional hazards: the
variance of the curves is not constant over time.
}
\examples{
data(bmt, package="KMsurv")
bmt <- within(bmt, {
z1 <- z1 -28
z2 <- z2- 28
z3 <- z1*z2
z4 <- as.double( group== 2 )
z5 <- as.double( group== 3 )
z6 <- z8
z7 <- (z7 / 30) - 9
z8 <- z10
})
c1 <- coxph(Surv(t2,d3) ~ z1+z2+z3+z4+z5+z6+z7+z8,
            method="breslow",
            data=bmt)
dxPlot(c1)
data(alloauto, package="KMsurv")
c1 <- coxph(Surv(time,delta) ~ factor(type),
            method='breslow',
            data=alloauto)
dxPlot(c1)
c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
dxPlot(c1)
}
\references{
Examples are from Klein J, Moeschberger M 2003
\emph{Survival Analysis}, 2nd edition. New York: Springer.
Example 11.1 - 11.7, pg 355-66

Last example is from: Therneau T, Grambsch P 2000.
\emph{Modeling Survival Data}, 1st edition. New York:
Springer. Section 5.1.2, pg 91.

Andersen PK, Borgan O, Gill R, Keiding N 1982. Linear
Nonparametric Tests for Comparison of Counting Processes,
with Applications to Censored Survival Data, Correspondent
Paper. \emph{International Statistical Review}
\bold{50}(3):219--44.
\href{http://www.jstor.org/stable/1402489}{JSTOR}
}

