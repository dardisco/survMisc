\name{gamTerms}
\alias{gamTerms}
\title{Individual terms of a Generalized Additive or Cox Proportional
Hazard Model}
\usage{
  gamTerms(fit, se.fit = TRUE, link = FALSE, weights)
}
\arguments{
  \item{fit}{The result of a Generalized Additive Model
  (\code{gam}).}

  \item{se.fit}{if \code{TRUE}, return the standard errors
  as well.}

  \item{link}{if \code{TRUE}, then the individual terms are
  centered so that the average of the inverse-link of the
  data, i.e., the data on the original scale, has mean
  zero.}

  \item{weights}{a vector of case weights, data is centered
  so that the weighted mean is zero}
}
\value{
  A list with one element per term. Each element is a
  matrix whose columns are x, y, and (optionally) se(y),
  with one row per unique x value, sorted by the first
  column. (This makes it easy to plot the results). The
  first element of the list, \code{constant}, contains an
  overall mean for the decomposition.
}
\description{
  Returns the individual terms of a \code{gam} or
  \code{coxph} object, along with the standard errors, in a
  way useful for plotting.
}
\examples{
data(air)
gfit <- gam::gam(ozone ~ s(temperature) + s(wind), data=air)
temp <- gamTerms(gfit)
identical( names(temp), c("constant", "temperature", "wind") )
dim(temp$wind)    # air has 111 rows, but only 28 unique wind speeds
### plot the fit versus square root of wind speed
yy <- cbind(temp$wind[,2], temp$wind[,2] - 1.96*temp$wind[,3], temp$wind[,2] + 1.96*temp$wind[,3])
### Adding the constant makes this a plot of actual y (ozone) at
### the mean temp
yy <- yy + temp$constant
graphics::matplot(sqrt(temp$wind[,1]), yy, lty=c(1,2,2),
type='l', col=1, xaxt='n', xlab='Wind Speed', ylab='Ozone')
temp <- seq(3,19,2)
graphics::axis(1, sqrt(temp), format(temp))
}
\author{
  Terry Therneau, Dirk Larson, updated from S-plus by Chris
  Dardis
}
\seealso{
  \code{gam}, \code{plot.gam}
}
\keyword{plot}

