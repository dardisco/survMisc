\name{mean.Surv}
\alias{mean.Surv}
\title{Mean for \code{Surv} object}
\usage{
  \method{mean}{Surv} (x, alpha = 0.05, method = "Efron",
    tEnd = "", by = 1, dfm = FALSE, ...)
}
\arguments{
  \item{x}{A \code{Surv} object}

  \item{alpha}{Significance level \eqn{\alpha}{alpha}}

  \item{method}{If the last observation is censored at time
  \eqn{t_k}{t_k}, one of the following values for
  \eqn{\hat{S}}{S}, the Kaplan-Meier estimate of survival
  time from then until \code{tEnd} is used: \describe{
  \item{Effron}{ \eqn{\hat{S}=0}{S=0}} \item{Gill}{
  \eqn{\hat{S}=\hat{S}(t_k)}{S = S(t_k)}. I.e.
  \eqn{\hat{S}}{S} is the equal to the last recorded value
  of \eqn{\hat{S}}{S}}. \item{Brown}{ \eqn{\hat{S}=
  e^{\frac{t_i}{t_k} \log{\hat{S}(t_k)}} }{ e^[ (t_i/t_k)
  log S(t_k) ]} for \eqn{ t_k \le t_i \le \code{tEnd}}{
  t(k) <= t(i) <= tEnd} }}}

  \item{tEnd}{If the last observation is censored at time
  \eqn{t_k}{t_k}, an estimate of \eqn{\hat{S}}{S} will be
  generated from \eqn{t_k}{t_k} to \code{tEnd}. \cr If
  \code{tEnd=""} a value of \eqn{2t_{max}}{2*tMax}, twice
  the longest time recorded, is used.}

  \item{by}{Increments (units of time) between
  \eqn{t_k}{t_k} and \code{tEnd}}

  \item{dfm}{If \code{TRUE}, will return the dataframe used
  to calculate the statistics}

  \item{...}{Additional arguments}
}
\value{
  A list with the following elements: \item{mean}{Mean of
  the \code{Surv} object} \item{variance}{The variance}
  \item{CI}{The confidence level (from \eqn{alpha}{alpha}
  above)} \item{upper}{Upper value for confidence interval}
  \item{lower}{Lower value for the confidence interval} If
  the last observation is censored at time \eqn{t_k}{t_k},
  two values are returned, one calculated up to
  \eqn{t_k}{t_k}, the other to \code{tEnd}.
}
\description{
  Mean for \code{Surv} object
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
mean(s1)
s1 <- Surv(time=c(6,14,21,44,62),event=c(1,1,0,1,1))
mean(s1)
s2 <- Surv(time=c(6,14,21,44,62),event=c(1,1,0,1,0))
mean(s1)
}

