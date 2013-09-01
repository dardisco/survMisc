\name{survEst}
\alias{Gw}
\alias{KM}
\alias{NelA}
\alias{NelAVar}
\alias{survEst}
\title{Estimates of survival functions based on n and e}
\usage{
  survEst(n, e)

  KM(n, e)

  Gw(n, e)

  NelA(n, e)

  NelAVar(n, e)
}
\arguments{
  \item{n}{Vector with corresponding no. at risk at time
  points}

  \item{e}{Vector with no. events (e.g. deaths) at time
  points}
}
\value{
  \code{KM} returns \eqn{\hat{S}}{S}, the Kaplan-Meier
  estimator \cr

  \code{Gw} returns \eqn{Var [\hat{S}(t)]}, the variance of
  the Kaplan-Meier estimator \cr

  \code{NelA} returns \eqn{\hat{H}(t)}, the Nelson-Aalen
  estimator \cr

  \code{NelAV} returns \eqn{Var \hat{H}(t)}, the variance
  of the Nelson-Aalen estimator
}
\description{
  Estimates of survival functions based on n and e
}
\details{
  \code{KM} gives the Kaplan-Meier (product-limit)
  estimator of survival function. Value is \eqn{1} at
  \eqn{t=0} and thereafter given by: \deqn{\prod_{t \leq
  t_i} (1-\frac{e_i}{n_i} )}{PROD 1-e/n}

  \code{Gw} gives Greenwoods estimtor of the Kaplan-Meier
  (product-limit) estimator of survival function is given
  by:

  \deqn{Var[\hat{S}(t)] = [\hat{S}(t)]^2 \sum_{t_i \leq t}
  \frac{ \frac{e_i}{n_i}}{n_i(1- \frac{e_i}{n_i}) } }{
  Var[S](t) = S(t)^2 SUM (e/n) / n(1-n) }

  \code{NelA} gives the Nelson Aalen estimator of hazard
  function. Value is \eqn{\hat{H}=0}{H=0} at \eqn{t=0} and
  thereafter given by: \deqn{\hat{H}(t) = \sum_{t \leq t_i}
  \frac{e_i}{n_i} }{ H(t) = SUM e/n}

  \code{NelAVar} gives the variance of the Nelson-Aalen
  estimator, given by: \deqn{ Var[\hat{H}(t))] = \sum_{t_i
  \leq t} \frac{e_i}{n_i^2} }{ SUM e/(n^2)}
}
\note{
  \code{survEst} is a generic name for the functions
  documented.
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
t1 <- tne(s1)
KM(n=t1$n,e=t1$e)
###
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
t1 <- tne(s1)
Gw(n=t1$n,e=t1$e)
###
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
t1 <- tne(s1)
NelA(e=t1$e,n=t1$n)
###
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
t1 <- tne(s1)
NelAVar(n=t1$n,e=t1$e)
}

