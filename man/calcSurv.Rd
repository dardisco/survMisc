\name{calcSurv}
\alias{calcSurv}
\title{Estimators of survival & hazard functions with variance for a \code{Surv} object describing right censored data}
\usage{
  calcSurv(s)
}
\arguments{
  \item{s}{A \code{Surv} object}
}
\value{
  A data frame. \cr \cr Rows - one for each time where an
  event occurs. \cr \cr Columns are: \item{t}{time}
  \item{n}{no. at risk} \item{e}{no. events}
  \item{KM}{Survival estimate by Kaplan-Meier
  (Product-Limit) estimator} \item{KMV}{Variance of
  Kaplan-Meier estimate (Greenwoods formula)} \item{SNel}{
  Survival estimate from Nelson-Aalen estimator:
  \eqn{\hat{S}=e^{\hat{H}} }{S=e^H} }
  \item{HNel}{Nelson-Aalen estimate of hazard function}
  \item{HNelV}{Variance of Nelson-Aalen estimate}
  \item{HKM}{ Hazard estimate from Kaplan-Meier estimator:
  \eqn{\hat{H}=-log{\hat{S} }}{H = -log(S)} }
}
\description{
  Estimators of survival & hazard functions with variance
  for a \code{Surv} object describing right censored data
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
calcSurv(s1)
}
\seealso{
  \code{\link{tne.Surv}}

  \code{\link{KM}}

  \code{\link{Gw}}

  \code{\link{NelA}}

  \code{\link{NelAVar}}
}

