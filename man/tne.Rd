\name{tne}
\alias{tne}
\alias{tne.Surv}
\alias{tne.survfit}
\title{Time, No. at risk, No. events}
\usage{
  tne(x, ...)

  \method{tne}{Surv} (x, ...)

  \method{tne}{survfit} (x, ...)
}
\arguments{
  \item{x}{A object of class \code{Surv} or \code{survfit}}

  \item{...}{Additional arguments, e.g. \cr onlyEvents (For
  \code{Surv} object): if \code{TRUE} shows only times at
  which at least one event occurred. \cr Otherwise shows
  \emph{all} times recorded (including those censored).}
}
\value{
  For a \code{Surv} object: A data frame with columns:
  \item{t}{time} \item{n}{no. at risk} \item{e}{no. events}
  \cr \cr For a \code{survfit} object: A data frame with
  columns: \item{t}{time} \item{n}{no. at risk (total)}
  \item{e}{no. events (total)} \item{p}{predictor}
  \item{np}{no. at risk (by predictor)} \item{Ep}{no.
  events expected (by predictor)} \item{e_Ep}{no. events
  minus no. events expected}
}
\description{
  For a \code{Surv} object, gives time, no. at risk and no.
  events. \cr \cr For a \code{survfit} object with one
  predictor object describing right censored data, gives
  time, no. at risk and no. events. Also per-predictor
  estimates of no. at risk, no. expected and no. events -
  no. expected. \cr \cr No. events expected (per predictor)
  is given by: \deqn{\frac{e_i(n[p]_i)}{n_i}}{ e(i)n[p](i)
  / n(i)} where \eqn{n[p]_i} is the no. at risk for the
  predictor.
}
\details{
  tne
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
s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
tne(s)
s <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
tne(s)
}

