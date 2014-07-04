\name{plotTerm}
\alias{plotTerm}
\title{Plot individual terms of a Generalized Additive \code{gam}
or Cox Proportional Hazards \code{coxph} Model}
\usage{
  plotTerm(x, term = 1, se = TRUE, p = 0.95, rug = TRUE,
    const = 0, col = 1, ...)
}
\arguments{
  \item{x}{The result of a Generalized Additive Model
  (\code{gam}) or a Cox proportional hazards model
  (\code{coxph})}

  \item{term}{The term to be plotted}

  \item{se}{if \code{se=TRUE}, also plot confidence
  intervals based on the standard errors}

  \item{p}{P-value used to plot the confidence intervals
  from standard errors. Based on normal distribution}

  \item{rug}{Add rug (1-dimensional plot) to x-axis. See
  \code{?graphics::rug}}

  \item{const}{Value for constant term. If
  \code{const=TRUE}, add the overall mean for the
  decomposition as returned by \code{\link{gamTerms}}}

  \item{col}{Color of line(s) on plot. If \code{se=TRUE},
  use a vector of 3 colors: first is main line, second is
  lower CI, third upper CI.}

  \item{...}{Additional arguments passed to
  \code{graphics}}
}
\value{
  A plot (base graphics) of the term in question. If
  \code{se=TRUE}, this is done using
  \code{graphics::matplot} otherwise \code{graphics::plot}
  is used.
}
\description{
  Plot individual terms of a Generalized Additive
  \code{gam} or Cox Proportional Hazards \code{coxph} Model
}
\examples{
fit1 <- coxph(Surv(time, status) ~ sex + pspline(age), data=lung)
plotTerm(fit1, term=2, rug=FALSE, ylab="Log-hazard",
 col=c("blue", "red", "red"))
}
\author{
  Terry Therneau. Updated from S-plus by Chris Dardis
}
\seealso{
  \code{\link{gamTerms}}
}

