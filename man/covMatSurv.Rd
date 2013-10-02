\name{covMatSurv}
\alias{covMatSurv}
\title{Covariance matrix for survival data}
\usage{
  covMatSurv(t, n, e, n1)
}
\arguments{
  \item{t}{time}

  \item{n}{no. at risk}

  \item{e}{no. events}

  \item{n1}{no. at risk (by predictor). \cr If 2 groups,
  should be given as a vector n1 with no. at risk for n1.
  \cr If more than 2 groups, a matrix with a column for
  each group.}
}
\value{
  An array. First two dimensions = no. groups. \cr Third
  dimension is no. observations (time points). \cr \cr
  Where there are two groups, the resulting sparse square
  matrix at time \eqn{i} has diagonal elements: \deqn{v_i =
  - \frac{n0_i n1_i e_i (n_i-e_i)}{n_i^2(n_i-1)}}{
  n0(i).n1(i).e(i).(n(i)-e(i)) / n(i)^2.(n(i) -1)} where
  \eqn{n1} is the no. at risk in group 1. \cr \cr For more
  than two groups, the resulting square matrix has diagonal
  elements: \deqn{ v_{kki} =
  \frac{n_{ki}(n_i-n_{ki})e_i(n_i-e_i)}{n_i^2(n_i-1)}}{
  v[k,k,i] = n[k](i).[n(i) - n[k](i)].e(i).[n(i) -e(i)] /
  n(i)^2.(n(i)-1)} and off diagonal elements: \deqn{
  v_{kli} = \frac{ -n_{ki}n_{li}
  e_i(n_i-e_i)}{n_i^2(n_i-1)}}{ v[k,l,i] = -
  n[k](i).n[l](i).e(i).[n(i)-e(i)] / n(i)^2.[n(i)-1]}
}
\description{
  Gives variance-covariance matrix for comparing survival
  data for two or more groups. Inputs are vectors
  corresponding to observations at a set of discrete time
  points for right censored data, except for \eqn{n1}, the
  no. at risk by predictor. This should be specified as a
  vector for one group, otherwise as a matrix with each
  column corresponding to a group.
}
\examples{
\dontrun{data(tneKidney)
covMatSurv(t=tneKidney$t,n=tneKidney$n,e=tneKidney$e,n1=tneKidney$n_1)
}
}
\seealso{
  Called by \code{\link{comp}}
}

