\name{genDfSurv}
\alias{genDfSurv}
\title{Generate data frame for survival analysis}
\usage{
  genDfSurv(f = 0, b = 2, c = 1, n = 100, pb = 0.5, nf = 3,
    rc = 0.8, pe = 0.5, t0 = 1, tmax = 100,
    asFactor = FALSE, model = TRUE, timelim = 5)
}
\arguments{
  \item{b}{\dfn{binomial predictors}, the number of
  predictors which are binary, i.e. limited to 0 or 1}

  \item{f}{\dfn{factors}, the number of predictors which
  are factors}

  \item{c}{\dfn{continuous predictors}, the number of
  predictors which are continuous}

  \item{n}{number of observations in the data frame}

  \item{nf}{the no. of levels in a factor}

  \item{pb}{\dfn{probability for binomnial predictors} the
  probability of binomial predictors being =1. E.g. if
  \code{pb=0.3}, 30\% will be 1s, 70\% will be 0s}

  \item{rc}{\dfn{ratio for continuous variables} the ratio
  of levels of continuous variables to the total number of
  observations \dfn{n} e.g. if \code{rc=0.8} and
  \code{n=100}, it will be in the range 1-80}

  \item{pe}{\dfn{probability of event} the probability of
  events (typically death/failure) occurring, i.e. are =1.
  E.g. if \code{pe=0.5}, 50\% will be 1s, 50\% will be 0s}

  \item{t0}{Lowest (starting) time}

  \item{tmax}{Highest (final) time}

  \item{asFactor}{if \code{TRUE}, predictors given as
  factors will be converted to factors in the data frame
  before the model is fit}

  \item{timelim}{function will timeout after \code{timelim}
  secs}

  \item{model}{If \code{TRUE} will return fitted model}
}
\value{
  A list with the following values: \item{df}{data frame
  with predictors (labelled \eqn{x1,x2, ..., xn}) and
  outcome (y), with \emph{n} rows (observations)}
  \item{model}{if \code{model = TRUE} model fit with \code{
  survival::coxph()} }
}
\description{
  Generates a data frame with a binary outcome (events),
  times and optional predictor variables. An optional coxph
  model is fitted to the data. Model is fitted with
  \code{survival::coxph()}
}
\note{
  Using \code{asFactor=TRUE} with factors which have a
  large number of levels (e.g. \code{nf >30}) on large
  datasets (e.g. \eqn{n >1000}) can cause fitting to be
  slow.
}
\examples{
set.seed(1)
genDfSurv()
genDfSurv(b=0, c=2, n=100, pe=0.7)
genDfSurv(b=1, c=0, n=1000)
genDfSurv(f=1, nf=4, b=1, c=0, asFactor=TRUE)
}
\keyword{datagen}

