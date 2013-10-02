\name{autoplot}
\alias{autoplot}
\alias{autoplot.survfit}
\title{Generate a ggplot for \code{survfit} object}
\usage{
  \method{autoplot}{survfit} (object, ..., alpha = 0.5,
    shape = 3, xlab = "Time", ylab = "Survival",
    title = "Marks show times with censoring", CI = FALSE,
    bands = FALSE, pval = FALSE, plotTable = FALSE,
    divideTime = 1, returnTable = FALSE)
}
\arguments{
  \item{object}{An object of class \code{survfit}}

  \item{...}{Additional arguments}

  \item{alpha}{Transparency to use for confidence intervals
  or bands}

  \item{shape}{Shape of marks to indicate censored
  onservations. \cr Default is 3 which gives vertical
  ticks. \cr Use 10 for circular marks.}

  \item{xlab}{Label for x axis}

  \item{ylab}{Label for y axis}

  \item{title}{Title for graph}

  \item{CI}{Include confidence intervals (plotted as
  lines). \cr These are taken from the \code{survfit}
  object.}

  \item{bands}{Include confidence bands (plotted as ribbons
  i.e. filled and joined by diagonal lines rather than
  steps)}

  \item{pval}{Add \eqn{p} value from log-rank test to lower
  left of plot}

  \item{plotTable}{Add table below plot. This contains one
  row for each time increment for each strata, giving the
  number at risk at that time.}

  \item{divideTime}{Divide time by this number to give time
  increments. A larger number may be preferable with large
  data sets to avoid crowding.}

  \item{returnTable}{Return \code{data.frame} used for
  table.}
}
\value{
  A \code{ggplot} and optionally a \code{data.frame} as
  above.
}
\description{
  Uses \code{ggplot2} to plot survival curves (Kaplan-Meier
  plot)
}
\examples{
data(kidney, package="KMsurv")
s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
autoplot(s1)
autoplot(s1, CI=TRUE, pval=TRUE, plotTable=TRUE, divideTime=5)
s1 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
autoplot(s1)
data(rectum.dat, package="km.ci")
s1 <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
### change confidence intervals to log Equal-Precision confidence bands
km.ci::km.ci(s1, method="logep")
autoplot(s1, bands=TRUE)
}
\author{
  Chris Dardis. Based on existing work by R. Saccilotto,
  Abhijit Dasgupta, Gil Tomas and Mark Cowley.
}

