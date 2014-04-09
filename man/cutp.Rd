\name{cutp}
\alias{cutp}
\alias{cutp.coxph}
\alias{cutp.survfit}
\title{Cutpoint for a continuous variable in a \code{coxph} or \code{survfit} model}
\usage{
cutp(x, ...)

\method{cutp}{coxph}(x, ..., var = "", plot = FALSE)

\method{cutp}{survfit}(x, ..., var = "", plot = FALSE)
}
\arguments{
  \item{x}{A \code{survfit} or \code{coxph} object}

  \item{...}{Additional arguments}

  \item{var}{Variable to test. Must be continuous (i.e. >2
  unique values)}

  \item{plot}{If \code{TRUE} will plot cut points against
  the statistic \eqn{S}}
}
\value{
A \code{data.frame} with columns: \item{CP}{The cut point.
Optimum value at which to divide the groups into those
\eqn{\geq}{>=} the cutpoint and those below} \item{Q}{The
test statistic} \item{p}{p-value} If \code{plot} =
\code{TRUE} a plot of cut points against values of the
log-rank test statistic \eqn{U}.
}
\description{
Cutpoint for a continuous variable in a \code{coxph} or
\code{survfit} model
}
\details{
The statistic is based on the score test from the Cox
model. For a cut point \eqn{\mu}{mu}, of a predictor vector
\eqn{G} the data is split into two groups, those \eqn{\geq
\mu}{>= mu} and those \eqn{< \mu}{< mu}. \cr \cr The
log-rank statistic \eqn{U} is calculated for each unique
element \eqn{k} in \eqn{G}: \deqn{U_k = \sum_{i=1}^D (
e_i^+ - n_i^+ \frac{e_i}{n_i} )}{ U_k = sum ( e1[i] -
n1[i].e[i]/n[i] ) } Where \eqn{e_i^+}{e1[i]} and
\eqn{n_i^+}{n1[i]} refer to the no. of events and no. at
risk in those above the cutpoint. \cr The sum is taken to
across distinct times with observed events, to \eqn{D}, the
largest of these. \cr It is normalized (standardized) in
the case of censoring by finding \eqn{\sigma^2}{s^2} which
is: \deqn{ \sigma^2 = \frac{1}{D-1} \sum_i^D ( 1 -
\sum_{j=1}^i \frac{1}{D+1-j} )^2 }{ s^2 = 1/(D-1) SUM[i to
D] { 1 - SUM[j to i] (1/(D-j+1))}^2 } The test statistic is
then \deqn{Q = \frac{\max |U_k|}{\sigma \sqrt{D-1}} }{ Q =
[ max |U_k| ] / [ s.(D-1)^0.5 ] } Under the null hypothesis
that the chosen cut-point does not predict survival, the
distribution of \eqn{Q} had a limiting distibution of the
supremum of the absolute value of a Brownian bridge, that
is \deqn{ p= Pr ( \sup Q \geq q ) = 2 \sum_{i=1}^{\infty}
(-1)^{i+1} \exp (-2 i^2 q^2) }{ P(Q >= q) = 2 SUM [i to
Inf] (-1)^(i+1).e^(-2.i^2.q^2) }
}
\examples{
data(kidtran, package="KMsurv")
k1 <- kidtran
k2 <- k1[k1$gender==1 & k1$race==2, ]
c1 <- coxph(Surv(time, delta) ~ age, data = k2)
cutp(c1, var="age", plot=TRUE)
k2 <- k1[k1$gender==2 & k1$race==2, ]
c1 <- coxph(Surv(time, delta) ~ age, data = k2)
cutp(c1, var="age")
}
\references{
Examples are from Klein J, Moeschberger M 2003
\emph{Survival Analysis}, 2nd edition. New York: Springer.
Example 8.3, pp 273-274.

Contal C, O'Quigley J, 1999 An application of changepoint
methods in studying the effect of age on survival in breast
cancer. \emph{Computational Statistics & Data Analysis}
\bold{30}(3):253--70.
\href{http://www.sciencedirect.com/science/article/pii/S0167947398000966}{ScienceDirect}
}

