\name{calcSurv}
\alias{calcSurv}
\title{Calculate estimators of survival & hazard functions
with variance for a \code{Surv} object describing right censored data}
\usage{
  calcSurv(s, round = 7)
}
\arguments{
  \item{s}{A \code{Surv} object}

  \item{round}{No. digits to which to round (for display)}
}
\value{
  A data frame. \cr \cr Rows - one for each time where an
  event occurs. \cr \cr Columns are: \item{t}{time}
  \item{n}{no. at risk} \item{e}{no. events}
  \item{KM}{Survival estimate by Kaplan-Meier
  (Product-Limit) estimator} \item{KMV}{Variance of
  Kaplan-Meier estimate (Greenwoods formula)} \item{SNelA}{
  Survival estimate from Nelson-Aalen estimator:
  \eqn{\hat{S}=e^{\hat{H}} }{S=e^H} }
  \item{HNel}{Nelson-Aalen estimate of hazard function}
  \item{HNelV}{Variance of Nelson-Aalen estimate}
  \item{HKM}{ Hazard estimate from Kaplan-Meier estimator:
  \eqn{\hat{H}=-log{\hat{S} }}{H = -log(S)} }
}
\description{
  Calculate estimators of survival & hazard functions with
  variance for a \code{Surv} object describing right
  censored data
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
calcSurv(s1)
data.frame(calcSurv(s1, round=4)[1:4],
          round(sqrt(calcSurv(s1)[5]),4),
          calcSurv(s1, round=4)[7],
          round(sqrt(calcSurv(s1)[8]),4)
          )
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==2, ] # AML patients
s1 <- Surv(time=b1$t2, event=b1$d3)
calcSurv(s1)
}
\references{
  First example is from: Klein J, Moeschberger M 2003
  \emph{Survival Analysis}, 2nd edition. New York:
  Springer. Table 4.3, pg 97.
}
\seealso{
  \code{\link{tne}}

  \code{\link{KM}}

  \code{\link{Gw}}

  \code{\link{NelA}}

  \code{\link{NelAVar}}
}

