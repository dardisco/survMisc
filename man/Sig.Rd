\name{Sig}
\alias{Sig}
\alias{Sig.coxph}
\title{Sig}
\usage{
  Sig(x, ...)

  \method{Sig}{coxph} (x, ...)
}
\arguments{
  \item{x}{A model of class \code{coxph}}

  \item{...}{Additional arguments}
}
\value{
  A \code{data.frame} with one fow for each coefficient in
  the original model. There are 3 columns, one for each of
  the tests: \itemize{ \item Wald \item LR (likelihood
  ratio) \item Score }
}
\description{
  Sig

  Significance tests of coefficients in survival model.
  These are: \describe{ \item{Wald test}{the statistic is:
  \deqn{ \frac{\hat{B}}{\hat{SE}}}{ Bhat/SE} where
  \eqn{\hat{B}}{Bhat} is the estimate of the coefficient
  and \eqn{\hat{SE}}{SE} is its standard error.  }
  \item{Likelihood Ratio test}{the statistic is the
  difference in the likelihood ratio of the original model
  and that with the coefficient omitted.  } \item{Score
  test}{Null hypothesis is that \eqn{\hat{B}=0}{Bhat=0}.
  The statistic is cacluated by refitting the model with
  the coefficient omitted to generate initial values. It is
  then fitted again with all covariates, using these values
  and setting \eqn{\hat{B}=0}{Bhat=0}.  } } All statistics
  are distributed as chi-square, with degrees of freedom =
  no. of coefficients \eqn{-1}.
}

