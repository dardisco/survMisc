\name{quantile}
\alias{median}
\alias{median.Surv}
\alias{median.coxph}
\alias{median.survfit}
\alias{quantile}
\alias{quantile.Surv}
\alias{quantile.coxph}
\alias{quantile.survfit}
\title{Quantiles and median for \code{Surv} and \code{survfit} objects}
\usage{
quantile(x, ...)

\method{quantile}{Surv}(x, ..., q = c(25, 50, 75), alpha = 0.05)

\method{quantile}{survfit}(x, ..., q = c(25, 50, 75), alpha = 0.05)

\method{quantile}{coxph}(x, ..., q = c(25, 50, 75), alpha = 0.05)

median(x, ...)

\method{median}{Surv}(x, ..., CI = FALSE, alpha = 0.05, method = c("log",
  "lin", "asr"))

\method{median}{survfit}(x, ..., CI = FALSE, alpha = 0.05,
  method = c("log", "lin", "asr"))

\method{median}{coxph}(x, ..., CI = FALSE, alpha = 0.05, method = c("log",
  "lin", "asr"))
}
\arguments{
  \item{x}{A \code{Surv} or \code{survfit} object}

  \item{...}{Additional arguments}

  \item{q}{(for \code{quantile}) Vector of quantiles
  (expressed as percentage).}

  \item{alpha}{Significance level \eqn{\alpha}{alpha}}

  \item{CI}{(for \code{median}) include confidence
  interval}

  \item{method}{(for \code{median}) method for confidence
  interval, as below. \cr One of \code{log}, \code{linear}
  or \code{arcsine-square root}.}
}
\value{
For \code{quantile}: \cr For a \code{Surv} object, a
\code{data.frame} with quantile, and upper and lower
confidence intervals using the 3 methods. For a
\code{survfit} object, a list with one element for each
stratum. Each element is a \code{data.frame} as above. \cr
\cr For \code{median}: \cr For a \code{Surv} object, the
median time. For a \code{survfit} object a
\code{data.frame} with one row for each stratum. \cr If
\code{CI} = \code{TRUE} then upper and lower confidence
intervals with one method (default is \code{log}).
}
\description{
Quantiles and median for \code{Surv} and \code{survfit}
objects
}
\details{
Confidence intervals are calculated from
\eqn{\sigma}{sigma} which is: \deqn{\sigma (t) = \sqrt{
\frac{Var[\hat{S}(t)]}{\hat{S}^2(t)}} }{ ( Var[S(t)] /
(S(t)^2) )^0.5} The intervals are: \cr

\itemize{

\item linear \deqn{ \hat{S}(t) \pm Z_{1- \alpha} \sigma (t)
\hat{S}(t)}{ S(t)+- Z(1-alpha) sigma(t) S(t)} Where
\eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival
estimate. \cr \cr

\item log transform \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}},
\hat{S}(t)^{\theta} ] }{ [S(t)^(1/theta), S(t)^theta]}
Where \eqn{\theta}{theta} is: \deqn{ \exp{ \frac{Z_{1-
\alpha} \sigma (t)}{ \log{\hat{S}(t)}}} }{ exp (
Z(1-alpha)sigma(t) / log(S(t)) )} \cr \cr

\item Arcsine-square root transform. \cr Upper: \deqn{
\sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}} - \frac{Z_{1-
\alpha}\sigma(t)}{2} \sqrt{
\frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{ sin^2( max[0, arcsin
S(t)^0.5 - Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])} Lower:
\deqn{ \sin^2(\min[\frac{\pi}{2},
\arcsin{\sqrt{\hat{S}(t)}} + \frac{Z_{1-
\alpha}\sigma(t)}{2} \sqrt{
\frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{ sin^2( min[pi/2,
arcsin S(t)^0.5 + Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}

}
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
quantile(s1)
b1 <- bmt[bmt$group==2, ] # AML low-risk patients
s1 <- Surv(time=b1$t2, event=b1$d3)
quantile(s1)
b1 <- bmt[bmt$group==3, ] # AML high-risk patients
s1 <- Surv(time=b1$t2, event=b1$d3)
quantile(s1)
###
s1 <- survfit(Surv(t2, d3)~ group, data=bmt)
quantile(s1)
###
c1 <- coxph(Surv(t2, d3)~ group, data=bmt)
quantile(c1)
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
median(s1)
median(s1, CI=TRUE)
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- survfit(Surv(t2, d3)~ group, data=bmt)
median(s1)
median(s1, CI=TRUE, method="asr")
c1 <- coxph(Surv(t2, d3) ~ group, data=bmt)
median(c1)
}
\references{
Examples for quantiles are from: Klein J, Moeschberger M
2003 \emph{Survival Analysis}, 2nd edition. New York:
Springer. Example 4.2, pg 121.
}

