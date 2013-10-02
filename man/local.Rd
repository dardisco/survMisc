\name{local}
\alias{local}
\alias{locLR}
\alias{locLR.coxph}
\alias{locScore}
\alias{locScore.coxph}
\alias{locWald}
\alias{locWald.coxph}
\title{Local tests for a model}
\usage{
  locScore(x, ...)

  \method{locScore}{coxph} (x, ..., all = FALSE, hypo = "",
    ties = "breslow")

  locLR(x, ...)

  \method{locLR}{coxph} (x, ..., all = FALSE, hypo = "",
    ties = "breslow")

  locWald(x, ...)

  \method{locWald}{coxph} (x, ..., all = FALSE, hypo = "")
}
\arguments{
  \item{x}{A model of class \code{coxph}}

  \item{...}{Additional arguments}

  \item{ties}{Method of handling ties when refitting model.
  Must be one of \code{breslow} or \code{efron}}

  \item{all}{Fit \emph{all} combinations of predictors}

  \item{hypo}{Hypothesis to test. This should be a vector
  of \eqn{0}s and \eqn{1}s of the same length as the no. of
  coefficients in the model (with at least one 0 and one
  1). Zeros indicate coefficients to exclude, ones indicate
  coefficients to keep. \cr}
}
\value{
  For \code{locScore} a \code{list} with the following
  items: \cr \item{coefficients}{coefficients from refitted
  model(s)} \item{test}{test for hypothesis} \cr \cr For
  \code{locLR} and \code{locWald}, a \code{data.frame}
  showing the hypothesis being tested and the results of
  the test.
}
\description{
  Local tests for a model
}
\details{
  The null hypothesis is that some of the coefficients in
  the model are zero (\eqn{H_0 : \hat{B}_i=0, \quad i \geq
  1}{H0: Bi=0}) vs. the alternative that at least one of
  them is nonzero. \cr \cr All of these tests are
  distributed as chi-square with degrees of freedom =
  number of excluded coefficients. \cr \cr For the
  \bold{score} test, the model is fitted again with the
  coefficients of interest excluded. \cr A value for the
  remaining coefficients is obtained. Then the complete
  model is fit again using these new values as initial
  values for those remaining coefficients and using zero as
  the initial value for the excluded coefficients. \cr
  Values for the excluded coefficients are generated
  without iteration. (I.e. the first values calculated,
  with no convergence towards maximum likelihood
  estimators). \cr The test is: \deqn{ \chi_{SC}^2 = U^T
  I^{-1} U}{ SC = U I^-1 U} where \eqn{U} is the score
  vector and \eqn{I^{-1}}{I^-1} is the covariance or
  inverse of the information matrix. (These are given by
  \code{colSums(survival::coxph.detail(x)$score)} and
  \code{x$var} respectively). \cr \cr For the
  \bold{likelihood ratio} test, the model is also refit
  with the coefficients of interest excluded. The
  likelihood ratios from the full model and those with
  coefficients excluded are used to construct the test:
  \deqn{ \chi^2_{LR} = 2(LR_{full} - LR_{excluded})}{ LR =
  2*( LR[full] - LR[excluded] )} \cr \cr The \bold{Wald}
  chi-squared statistic is given by: \deqn{ \chi^2_W =
  \hat{B}^T I^{-1} \hat{B} }{ W = B I^-1 B } Where
  \eqn{\hat{B}}{B} is the vector of fitted coefficients
  (from the complete model) thought to be \eqn{=0}. \cr
  \eqn{I^{-1}}{I^-1} is composed of the corresponding
  elements from the covariance matrix of the model.
}
\examples{
data(larynx, package="KMsurv")
c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx)
locScore(c1, all=TRUE)
locScore(c1, hypo=c(0,0,0,1))
###
data(larynx, package="KMsurv")
c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx, method="breslow")
locLR(c1, all=TRUE)
locLR(c1, hypo=c(0,0,0,1))
###
data(larynx, package="KMsurv")
c1 <- coxph(Surv(time, delta) ~ factor(stage) + age, data=larynx,
method="breslow")
locWald(c1, all=TRUE)
locWald(c1, hypo=c(0,0,0,1))
}
\references{
  Examples are from: Klein J, Moeschberger M 2003
  \emph{Survival Analysis}, 2nd edition. New York:
  Springer. Example 8.2, pp 264-6.
}

