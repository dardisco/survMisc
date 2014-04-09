\name{tne}
\alias{tne}
\alias{tne.Surv}
\alias{tne.coxph}
\alias{tne.survfit}
\title{Time, No. at risk, No. events}
\usage{
tne(x, ...)

\method{tne}{Surv}(x, ..., onlyEvents = FALSE)

\method{tne}{survfit}(x, ..., onlyEvents = FALSE)

\method{tne}{coxph}(x, ..., onlyEvents = FALSE)
}
\arguments{
  \item{x}{A object of class \code{Surv} or \code{survfit}}

  \item{...}{Additional arguments, e.g. for \code{survfit}
  object: \cr \cr \code{pred} : if \code{TRUE} gives
  predicted values as below \cr \cr \code{asList} : if
  \code{TRUE} will return a list with one element for each
  stratum}

  \item{onlyEvents}{if \code{TRUE} shows only times at
  which at least one event occurred. Otherwise shows
  \emph{all} times recorded (including those censored).}
}
\value{
For a \code{Surv} object: A data frame with columns:
\item{t}{time} \item{n}{no. at risk} \item{e}{no. events}
\cr \cr For a \code{survfit} or \code{coxph} object: A data
frame with columns: \item{t}{time} \item{n}{no. at risk (by
strata)} \item{e}{no. events (by strata)} \item{s}{strata}
\cr \cr If \code{asList} = \code{TRUE} then instead a
\code{list} with one element for each stratum, where each
elements is a \code{data.frame} with columns \strong{t},
\strong{n} and \strong{e} as above. \cr \cr If \code{pred}
= \code{TRUE} then instead a \code{data.frame} with
columns: \item{t}{time} \item{n}{no. at risk (total)}
\item{e}{no. events (total)} \item{s}{strata} \item{ns}{no.
at risk (by strata)} \item{Es}{no. events expected (by
strata)} \item{e_Es}{no. events minus no. events expected}
}
\description{
For a \code{Surv} object, gives time, no. at risk and no.
events. \cr \cr For a \code{survfit} or \code{coxph} object
describing right censored data, gives time, no. at risk and
no. events by strata. \cr Also (optionally) per-stratum
estimates of no. at risk, no. expected and no. events - no.
expected. \cr \cr No. events expected (per predictor) is
given by: \deqn{\frac{e_i(n[p]_i)}{n_i}}{ e(i)n[p](i) /
n(i)} where \eqn{n[p]_i} is the no. at risk for the
predictor.
}
\examples{
### Surv object
df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
                  e=rep(c(0,1),4))
s1 <- Surv(df0$t,df0$e,type="right")
tne(s1)
tne(s1, onlyEvents=FALSE)
### survfit object
data(kidney, package="KMsurv")
s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
tne(s1)
tne(s1, asList=TRUE)
tne(s1, pred=TRUE)
s2 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
tne(s2)
data(larynx, package="KMsurv")
s1 <- survfit(Surv(time, delta) ~ factor(stage) + age, data=larynx)
tne(s1, pred=TRUE)
data(bmt, package="KMsurv")
s1 <- survfit(Surv(t2, d3) ~ z3 +z10, data=bmt)
tne(s1)
s1 <- survfit(Surv(t2, d3) ~ 1, data=bmt)
tne(s1, pred=TRUE, onlyEvents=TRUE)
### coxph object
data(kidney, package="KMsurv")
c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney)
tne(c1)
tne(c1, asList=TRUE)
}
\references{
Example using \code{kidney} data is from: Klein J,
Moeschberger M 2003 \emph{Survival Analysis}, 2nd edition.
New York: Springer. Example 7.2, pg 210.
}

