\name{sig}
\alias{sig}
\alias{Sig}
\alias{sig.coxph}
\title{Significiance tests of coefficients in a Coxph model}
\usage{
  Sig(x, ...)

  sig(x, ...)

  \method{sig}{coxph} (x, ...)
}
\arguments{
  \item{x}{A model of class \code{coxph}}

  \item{...}{Additional arguments (not implemented)}
}
\value{
  A \code{data.frame} with one fow for each coefficient in
  the original model. There are 3 columns, one for each of
  the tests: \itemize{ \item Wald \item LR (likelihood
  ratio) \item Score }
}
\description{
  These are: \describe{ \item{Wald test}{the statistic is:
  \deqn{ \frac{\hat{B}}{\hat{SE}}}{ Bhat/SE} where
  \eqn{\hat{B}}{Bhat} is the estimate of the coefficient
  and \eqn{\hat{SE}}{SE} is its standard error.  }
  \item{Likelihood Ratio test}{the statistic is the
  difference in the likelihood ratio of the original model
  and that with the coefficient omitted.  } \item{Score
  test}{Aka the \bold{log-rank} test. Null hypothesis is
  that \eqn{\hat{B}=0}{Bhat=0}. The statistic is cacluated
  by refitting the model with the coefficient omitted to
  generate initial values. It is then fitted again with all
  covariates, using these values and setting
  \eqn{\hat{B}=0}{Bhat=0}.  } } All statistics are
  distributed as chi-square, with degrees of freedom
  \eqn{=} no. of coefficients \eqn{-1}.
}

