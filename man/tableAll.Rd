\name{tableAll}
\alias{tableAll}
\title{Table all predictors in a formula with the outcome}
\usage{
tableAll(formula = y ~ ., data = parent.frame(), return = c("summary",
  "zeros", "zEq", "counts", "all"), nlf = 2)
}
\arguments{
  \item{formula}{A formula. Also works with formulas where
  the left-hand side is a \code{Surv} object describing
  right-censored data}

  \item{data}{A \code{data.frame}}

  \item{return}{See \bold{Value} below}

  \item{nlf}{Number of levels defining a factor. Predictors
  with \eqn{>nlf} levels are considered continuous and are
  not tabulated. Needs to be less than the number of rows
  in the \code{data.frame}}
}
\value{
\itemize{ \item If \code{return="summary"} (the default), a
table with one row per predictor and three columns:
\describe{ \item{zeros}{at least one zero present}
\item{someEq}{outcomes equal for least \emph{some} levels
of the predictor} \item{allEq}{outcomes equal for
\emph{all} levels of the predictor} } } Other values return
a list of tables. Each element is named after the
predictor. \itemize{ \item If \code{return="zeros"}, one
table for each predictor with a least one zero present. \cr
Each table shows only those levels of the predictor for
which one level of the outcome is zero.  \item If
\code{return="zEq"}, one table for each predictor with a
least one zero present or one level which has equal
outcomes. \cr Each table shows only those levels where one
of the above apply.  \item If \code{return="counts"}, each
tables gives the total number of levels where zeros and
equal outcomes are present and absent.  \item If
\code{return="all"}, a list of tables of outcomes for
\emph{all} levels of each predictor. }
}
\description{
Table all predictors in a formula with the outcome
}
\details{
Cross-tabulation of outcomes against levels of a predictor
is useful in models where the outcome has limited values.
}
\examples{
set.seed(1)
d1 <- genDfSurv(c=3, rc=0.5)$df
tableAll(Surv(t1, e) ~ ., data=d1, return="summary", nlf=99)
t1 <- tableAll(Surv(t1, e) ~ ., data=d1, return="c", nlf=99)
### simple graph
p <- par()
par( mfrow=c(2, 2))
for (i in 1:length(t1)){
    graphics::mosaicplot(t1[[i]], main="", cex=1.5)
    title(main=list(names(t1[i]), cex=3))
}
par <- p
set.seed(2)
d1 <- genDfSurv(f=1, n=10)$df
t1 <- tableAll(Surv(t1, e) ~ x1, nlf=9, data=d1)
tableAll(e ~ x1, nlf=9, r="zEq", data=d1)
tableAll(e ~ ., nlf=3, r="c", data=d1)
}
\author{
Chris Dardis
}

