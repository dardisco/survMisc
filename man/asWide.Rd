% Generated by roxygen2 (4.0.1): do not edit by hand
\name{asWide}
\alias{asLong}
\alias{asLong.ten}
\alias{asWide}
\alias{asWide.ten}
\title{Convert an object to "wide" or "long" form.}
\usage{
asWide(x, ...)

\method{asWide}{ten}(x, ...)

asLong(x, ...)

\method{asLong}{ten}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ten} or \code{pred}.}

\item{...}{Additional arguments (not implemented).}
}
\value{
A new \code{data.table} is returned,
with the data in 'wide' or 'long' format.
 \cr
There is one row for each time point.
 \cr
For a \code{ten} object generated from a \code{numeric} or \code{Surv} object,
 this has columns:
 \item{t}{\bold{t}ime.}
 \item{e}{number of \bold{e}vents.}
 \item{n}{\bold{n}umber at risk.}
If derived from a \code{survfit}, \code{coxph} or \code{formula} object,
there are additional columns for \code{e} and \code{n}
for \emph{each} covariate group.
}
\description{
Convert an object to "wide" or "long" form.






}
\note{
Most methods for \code{ten} objects are designed for the 'long' form.
}
\examples{
data("bmt", package="KMsurv")
require("survival")
t1 <- ten(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
asWide(t1)
asLong(asWide(t1))
stopifnot(asLong(asWide(t1)) == ten(ten(t1)))
}

